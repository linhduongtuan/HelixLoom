//! The PhyloCNN networks on candle (a Rust tensor library): the same layers as `Network::run`, many trees
//! at a time, on the GPU where there is one (Metal on a Mac, CUDA when built for it) or on candle's own CPU
//! code. The plain-Rust `Network::run` stays the reference; the tests compare the two.

use super::{Act, Layer, Network};
use candle_core::{Device, Tensor};

/// Trees put through the network together; bounds the memory the GPU needs.
const CHUNK: usize = 128;

/// The best GPU there is, with a name for the window; `None` when there is only the CPU.
pub fn device() -> Option<(Device, String)> {
    #[cfg(target_os = "macos")]
    if let Ok(d) = Device::new_metal(0) {
        return Some((d, "Apple GPU (Metal)".into()));
    }
    #[cfg(feature = "cuda")]
    if let Ok(d) = Device::new_cuda(0) {
        return Some((d, "NVIDIA GPU (CUDA)".into()));
    }
    None
}

fn bad(e: candle_core::Error) -> String {
    format!("GPU: {e}")
}

enum Step {
    Conv {
        kernel: Tensor,
        groups: usize,
        act: Act,
    },
    Norm {
        scale: Tensor,
        shift: Tensor,
    },
    Pool,
    Dense {
        w: Tensor,
        b: Tensor,
        act: Act,
    },
}

fn activate(x: &Tensor, act: Act) -> candle_core::Result<Tensor> {
    match act {
        Act::Elu => x.elu(1.0),
        Act::Relu => x.relu(),
        Act::Linear => Ok(x.clone()),
        Act::Softmax => {
            let e = x.broadcast_sub(&x.max_keepdim(1)?)?.exp()?;
            e.broadcast_div(&e.sum_keepdim(1)?)
        }
    }
}

/// The network's weights as tensors on `device`: convolution kernels as (filters, channels per group, 1,
/// width), batch normalisation as one scale and one shift per channel.
fn steps(net: &Network, device: &Device) -> candle_core::Result<Vec<Step>> {
    net.layers
        .iter()
        .map(|layer| {
            Ok(match layer {
                Layer::Conv {
                    groups,
                    kw,
                    cin,
                    filters,
                    act,
                    w,
                } => {
                    let per_group = cin / groups;
                    let mut kernel = vec![0.0f32; filters * per_group * kw];
                    for f in 0..*filters {
                        for c in 0..per_group {
                            for j in 0..*kw {
                                kernel[(f * per_group + c) * kw + j] =
                                    w[(j * per_group + c) * filters + f];
                            }
                        }
                    }
                    Step::Conv {
                        kernel: Tensor::from_vec(kernel, (*filters, per_group, 1, *kw), device)?,
                        groups: *groups,
                        act: *act,
                    }
                }
                Layer::Norm {
                    eps,
                    gamma,
                    beta,
                    mean,
                    var,
                } => {
                    let scale: Vec<f32> = gamma
                        .iter()
                        .zip(var)
                        .map(|(g, v)| g / (v + eps).sqrt())
                        .collect();
                    let shift: Vec<f32> = beta
                        .iter()
                        .zip(mean)
                        .zip(&scale)
                        .map(|((b, m), s)| b - m * s)
                        .collect();
                    let n = scale.len();
                    Step::Norm {
                        scale: Tensor::from_vec(scale, (1, n, 1, 1), device)?,
                        shift: Tensor::from_vec(shift, (1, n, 1, 1), device)?,
                    }
                }
                Layer::Pool => Step::Pool,
                Layer::Dense {
                    inputs,
                    outputs,
                    act,
                    w,
                    b,
                } => Step::Dense {
                    w: Tensor::from_vec(w.clone(), (*inputs, *outputs), device)?,
                    b: Tensor::from_vec(b.clone(), (1, *outputs), device)?,
                    act: *act,
                },
            })
        })
        .collect()
}

/// Put every input through the network on `device`; one output vector per input.
pub fn run_batch_on(
    net: &Network,
    inputs: &[Vec<f32>],
    device: &Device,
) -> Result<Vec<Vec<f32>>, String> {
    let [rows, width, channels] = net.input;
    if let Some(wrong) = inputs.iter().find(|x| x.len() != rows * width * channels) {
        return Err(format!(
            "The network takes {rows} x {width} x {channels} numbers, not {}.",
            wrong.len()
        ));
    }
    let layers = steps(net, device).map_err(bad)?;
    let mut out = Vec::with_capacity(inputs.len());
    for chunk in inputs.chunks(CHUNK) {
        let flat: Vec<f32> = chunk.iter().flatten().copied().collect();
        // the network's rows x width x channels become channels x rows x width, as candle's convolutions want
        let mut x = Tensor::from_vec(flat, (chunk.len(), rows, width, channels), device)
            .and_then(|t| t.permute((0, 3, 1, 2))?.contiguous())
            .map_err(bad)?;
        for layer in &layers {
            x = match layer {
                Step::Conv {
                    kernel,
                    groups,
                    act,
                } => x
                    .conv2d(kernel, 0, 1, 1, *groups)
                    .and_then(|y| activate(&y, *act)),
                Step::Norm { scale, shift } => {
                    x.broadcast_mul(scale).and_then(|y| y.broadcast_add(shift))
                }
                Step::Pool => x.mean((2, 3)),
                Step::Dense { w, b, act } => x
                    .matmul(w)
                    .and_then(|y| y.broadcast_add(b))
                    .and_then(|y| activate(&y, *act)),
            }
            .map_err(bad)?;
        }
        out.extend(x.to_vec2::<f32>().map_err(bad)?);
    }
    Ok(out)
}

/// The same on the best GPU; an error when there is none.
pub fn run_batch(net: &Network, inputs: &[Vec<f32>]) -> Result<Vec<Vec<f32>>, String> {
    let (device, _) = device().ok_or("This computer has no GPU the program can use.")?;
    run_batch_on(net, inputs, &device)
}
