"""HelixLoom's PyTorch worker: one process, one device (a CUDA card, the Apple GPU, or the CPU for tests).

The Rust program starts one of these for every GPU a job was given (`src/accel/torch/mod.rs`), sends it the columns of an alignment
once, and then asks, over its standard input and output, for the counts that make a matrix of distances and for the likelihood of
a tree. The program does the searching, the queueing and the summing; this does only the tensor work, in 32 bits, so that every
answer can be compared with the processor's.

Messages in both directions are binary: an operation name, then named arrays.

    u8 length, the name | u16 number of arrays | for each: u8 length, name | u8 type | u64 count | the bytes (little-endian)

Types: 1 u8, 2 u32, 3 f32, 4 f64. A reply is `ok` with its arrays, or `error` with a text. The standard output carries nothing else:
whatever PyTorch says goes to the standard error, which the job's log keeps.
"""

from __future__ import annotations

import os
import struct
import sys
import time
import traceback
from typing import Callable

import numpy as np
import torch

U8, U32, F32, F64 = 1, 2, 3, 4
KINDS = {U8: np.uint8, U32: np.dtype("<u4"), F32: np.dtype("<f4"), F64: np.dtype("<f8")}
CODES = {np.dtype(np.uint8): U8, np.dtype("<u4"): U32, np.dtype("<f4"): F32, np.dtype("<f8"): F64}

# What a tip says, a row for each of A, C, G, T and "any base": which states at the tip fit.
TIP_STATES = ((1.0, 0.0, 0.0, 0.0), (0.0, 1.0, 0.0, 0.0), (0.0, 0.0, 1.0, 0.0), (0.0, 0.0, 0.0, 1.0), (1.0, 1.0, 1.0, 1.0))

stdin = sys.stdin.buffer
# The protocol owns the program's standard output and nothing else may write there: it is a copy of the descriptor (a buffered writer,
# which writes all of a long reply; the unbuffered one a `-u` run gives can stop short), and descriptor 1 itself is now the standard
# error, so that anything a library prints, or C code writes, lands in the job's log and not in the middle of a reply.
stdout = os.fdopen(os.dup(1), "wb")
os.dup2(2, 1)
sys.stdout = sys.stderr


def read_exactly(n: int) -> bytes:
    data = stdin.read(n)
    if len(data) != n:
        raise EOFError
    return data


def read_message() -> tuple[str, dict[str, np.ndarray]]:
    """The next message as (operation, {name: array}); EOFError when the program has closed our input."""
    first = stdin.read(1)
    if not first:
        raise EOFError
    op = read_exactly(first[0]).decode()
    (count,) = struct.unpack("<H", read_exactly(2))
    arrays = {}
    for _ in range(count):
        name = read_exactly(read_exactly(1)[0]).decode()
        kind = read_exactly(1)[0]
        (n,) = struct.unpack("<Q", read_exactly(8))
        dtype = np.dtype(KINDS[kind])
        arrays[name] = np.frombuffer(read_exactly(n * dtype.itemsize), dtype=dtype)
    return op, arrays


def write_message(op: str, arrays: dict | None = None) -> None:
    parts = [bytes([len(op)]), op.encode()]
    arrays = arrays or {}
    parts.append(struct.pack("<H", len(arrays)))
    for name, value in arrays.items():
        if isinstance(value, str):
            value = np.frombuffer(value.encode(), dtype=np.uint8)
        value = np.ascontiguousarray(value)
        if value.dtype not in CODES:
            raise TypeError(f"cannot send {value.dtype}")
        parts += [bytes([len(name)]), name.encode(), bytes([CODES[value.dtype]]), struct.pack("<Q", value.size), value.tobytes()]
    stdout.write(b"".join(parts))
    stdout.flush()


def scalar(arrays: dict[str, np.ndarray], name: str) -> float:
    return arrays[name][0].item()


def pick_device(spec: str) -> tuple[torch.device, str, int]:
    """(device, kind, how many of that kind there are). `auto` takes the first CUDA card, else the Apple GPU, else nothing."""
    if spec == "auto":
        if torch.cuda.is_available():
            return torch.device("cuda", 0), "cuda", torch.cuda.device_count()
        if torch.backends.mps.is_available():
            return torch.device("mps"), "mps", 1
        return torch.device("cpu"), "none", 0
    device = torch.device(spec)
    if device.type == "cuda":
        return device, "cuda", torch.cuda.device_count()
    return device, device.type, 1


def synchronize(device: torch.device) -> None:
    if device.type == "cuda":
        torch.cuda.synchronize(device)
    elif device.type == "mps":
        torch.mps.synchronize()


def describe(device: torch.device, kind: str) -> tuple[str, int]:
    if kind == "cuda":
        p = torch.cuda.get_device_properties(device)
        return f"{p.name} ({p.total_memory >> 20} MB, sm_{p.major}{p.minor})", p.total_memory
    if kind == "mps":
        return "Apple GPU (Metal)", 0
    return kind.upper(), 0


class Worker:
    def __init__(self, device: torch.device, kind: str) -> None:
        self.device = device
        self.kind = kind
        self.sets = {}
        # the branches being scanned: what lies beyond both ends of each, for the patterns this worker was given
        self.edges = {}
        # no rounding of the inputs of a product to a short mantissa (TF32): the counts must stay whole numbers, and the likelihood exact to 32 bits
        torch.backends.cuda.matmul.allow_tf32 = False
        torch.backends.cudnn.allow_tf32 = False
        torch.set_float32_matmul_precision("highest")

    # ---- speed -----------------------------------------------------------------------------------------------------------------

    def rate(self, work: Callable[[], object]) -> float:
        """How many times a second `work` runs: the median of three runs of at least `SPAN` seconds each, after a warm-up of as long (a
        card's clocks take a while to rise, and a run of a millisecond measures the launching, not the card)."""
        span = 0.15

        def run() -> float:
            start = time.perf_counter()
            calls = 0
            while True:
                work()
                calls += 1
                synchronize(self.device)
                elapsed = time.perf_counter() - start
                if elapsed >= span:
                    return calls / elapsed

        run()
        return sorted(run() for _ in range(3))[1]

    def speed(self, _: dict[str, np.ndarray]) -> dict[str, np.ndarray]:
        """How fast this device runs the two kinds of work, in runs a second. The likelihood is limited by memory, the distances by
        matrix products, and two cards that are alike in one can differ in the other."""
        d = self.device
        v = torch.ones(1 << 20, 4, device=d)
        p = torch.ones(4, 4, device=d)
        a = torch.ones(2048, 2048, device=d)

        def likelihood() -> None:
            mixed = v @ p.T
            top = mixed.amax(1, keepdim=True)
            (mixed / top).log()

        return {"memory": np.array([self.rate(likelihood)]), "matmul": np.array([self.rate(lambda: a @ a)])}

    # ---- distances -------------------------------------------------------------------------------------------------------------

    def load_distances(self, a: dict[str, np.ndarray]) -> dict:
        """Write the rows (u8, `n` by `columns`) out as 0 and 1 for the three products, in chunks of columns, and keep them here."""
        n, columns = int(scalar(a, "n")), int(scalar(a, "columns"))
        rows = torch.from_numpy(a["rows"].reshape(n, columns).copy())
        lut = torch.full((256,), -1, dtype=torch.int64)
        for k, letter in enumerate("ACGT"):
            lut[ord(letter)] = k
            lut[ord(letter.lower())] = k
        lut = lut.to(self.device)
        chunk = max(64, min(8192, (256 << 20) // (n * 7 * 4)))
        chunks = []
        for start in range(0, columns, chunk):
            # the letters go to the card a chunk at a time and are written out as 0 and 1 there
            codes = lut[rows[:, start : start + chunk].to(self.device).long()]
            valid = (codes >= 0).float()
            base = codes.clamp(min=0)
            onehot = torch.nn.functional.one_hot(base, 4).float() * valid[..., None]
            purine = torch.nn.functional.one_hot(base % 2, 2).float() * valid[..., None]
            chunks.append((start, onehot, purine, valid))
        self.sets[int(scalar(a, "id"))] = ("distances", n, columns, chunks)
        # an error of the card shows here, not as the fault of the next request
        synchronize(self.device)
        return {}

    def counts(self, a: dict[str, np.ndarray]) -> dict[str, np.ndarray]:
        """For the rows `first` to `first + length` against all rows, the three products, each column taken `weights[column]` times."""
        _, n, columns, chunks = self.sets[int(scalar(a, "id"))]
        first, length = int(scalar(a, "first")), int(scalar(a, "length"))
        weights = torch.from_numpy(a["weights"].astype(np.float32)).to(self.device) if "weights" in a else None
        valid_both = torch.zeros(length, n, device=self.device)
        identical = torch.zeros_like(valid_both)
        same_class = torch.zeros_like(valid_both)
        for start, onehot, purine, valid in chunks:
            w = onehot.shape[1]
            ws = weights[start : start + w] if weights is not None else None
            for target, x in ((identical, onehot), (same_class, purine), (valid_both, valid)):
                x = x.reshape(n, -1) if x.dim() == 3 else x
                left = x[first : first + length]
                if ws is not None:
                    per = x.shape[1] // w
                    left = left * ws.repeat_interleave(per)[None, :]
                target += left @ x.T
        return {"valid": valid_both.cpu().numpy(), "identical": identical.cpu().numpy(), "same_class": same_class.cpu().numpy()}

    # ---- likelihood ------------------------------------------------------------------------------------------------------------

    def load_likelihood(self, a: dict[str, np.ndarray]) -> dict:
        tips, patterns = int(scalar(a, "tips")), int(scalar(a, "patterns"))
        states = torch.from_numpy(a["states"].reshape(tips, patterns).astype(np.int32)).to(self.device)
        counts = np.array(a["counts"], dtype=np.float64)
        self.sets[int(scalar(a, "id"))] = ("likelihood", tips, patterns, states, counts)
        synchronize(self.device)
        return {}

    def prune(self, steps: np.ndarray, inner_p: np.ndarray, tip_p: np.ndarray, states: torch.Tensor, first: int, length: int, first_node: int):
        """Felsenstein's pruning of the steps over the patterns `first` to `first + length`: the vector (patterns by 4) of the node
        `first_node`, scaled so that its largest number is 1, and the logarithms of the scales (patterns by 1)."""
        d = self.device
        steps = steps.reshape(-1, 7)
        # one matrix (or table) for each branch, taken apart once and not by an operation at every node
        inner = torch.from_numpy(inner_p.copy()).to(d).reshape(-1, 4, 4).unbind(0)
        tips = torch.from_numpy(tip_p.copy()).to(d).reshape(-1, 5, 4).unbind(0)
        done = {}
        for node, a_tip, a_node, a_slot, b_tip, b_node, b_slot in steps.tolist():
            product, scale = None, None
            for is_tip, child, slot in ((a_tip, a_node, a_slot), (b_tip, b_node, b_slot)):
                if is_tip:
                    said = torch.index_select(tips[slot], 0, states[child, first : first + length])
                else:
                    vector, below = done.pop(child)
                    said = vector @ inner[slot].T
                    scale = below if scale is None else scale + below
                product = said if product is None else product * said
            top = product.amax(1, keepdim=True)
            logs = top.log()
            done[node] = (product / top, logs if scale is None else scale + logs)
        return done.pop(first_node)

    def likelihood(self, a: dict[str, np.ndarray]) -> dict[str, np.ndarray]:
        """The sum over the patterns `first` to `first + length` of their counts times the logarithm of their likelihoods (in 64 bits
        here, from the 32-bit logarithms), by Felsenstein's pruning over all the patterns at once, as `src/seq/likelihood/gpu.rs` does."""
        _, tips, patterns, states, counts = self.sets[int(scalar(a, "id"))]
        first, length = int(scalar(a, "first")), int(scalar(a, "length"))
        d = self.device
        inner_p, tip_p = a["inner_p"], a["tip_p"]
        vector, scale = self.prune(a["steps"], inner_p, tip_p, states, first, length, int(scalar(a, "first_node")))
        root_p = torch.from_numpy(a["root_p"].copy()).to(d).reshape(4, 4)
        root_table = torch.from_numpy(a["root_table"].copy()).to(d).reshape(5, 4)
        said = vector @ root_p.T
        weights = torch.index_select(root_table, 0, states[0, first : first + length])
        logs = ((weights * said).sum(1).log() + scale.squeeze(1)).cpu().numpy().astype(np.float64)
        return {"sum": np.array([float(np.dot(counts[first : first + length], logs))])}

    # ---- the length of one branch ----------------------------------------------------------------------------------------------

    def edge_begin(self, a: dict[str, np.ndarray]) -> dict:
        """Keep, for the patterns `first` to `first + length` of a set, what lies beyond both ends of one branch: for each pattern, the
        likelihood of the part of the tree beyond one end given each state at that end, the base frequency of the state at the end `a`
        multiplied in, as the 16 products with the states at the end `b`. A branch's length is then tried by one product with its 16
        transition probabilities (`edge_eval`), not by the pruning of the whole tree, which is what the search does about twenty times for
        each branch it sets. A tip says what it says by a table of five rows; an inner end is the pruning of the steps given for it."""
        _, _, _, states, counts = self.sets[int(scalar(a, "id"))]
        first, length = int(scalar(a, "first")), int(scalar(a, "length"))
        d = self.device
        vectors, scales = [], []
        for end in ("a", "b"):
            node = int(scalar(a, end + "_node"))
            if int(scalar(a, end + "_tip")):
                table = torch.tensor(TIP_STATES, device=d)
                vectors.append(torch.index_select(table, 0, states[node, first : first + length]))
            else:
                vector, scale = self.prune(a[end + "_steps"], a[end + "_inner_p"], a[end + "_tip_p"], states, first, length, node)
                vectors.append(vector)
                scales.append(scale.squeeze(1))
        freqs = torch.from_numpy(a["freqs"].copy()).to(d)
        weights = (freqs[None, :, None] * vectors[0][:, :, None] * vectors[1][:, None, :]).reshape(-1, 16)
        scale = None
        for each in scales:
            scale = each if scale is None else scale + each
        self.edges[int(scalar(a, "scan"))] = (weights, scale, np.array(counts[first : first + length], dtype=np.float64))
        # an error of the card shows here, not as the fault of the next request
        synchronize(self.device)
        return {}

    def edge_eval(self, a: dict[str, np.ndarray]) -> dict[str, np.ndarray]:
        """The sum over the kept patterns of their counts times the logarithm of their likelihoods with the branch `p` (16 numbers: the
        matrix of its length) long, in 64 bits from the 32-bit logarithms, as `likelihood` does."""
        weights, scale, counts = self.edges[int(scalar(a, "scan"))]
        p = torch.from_numpy(a["p"].copy()).to(self.device)
        logs = (weights @ p).log()
        if scale is not None:
            logs = logs + scale
        return {"sum": np.array([float(np.dot(counts, logs.cpu().numpy().astype(np.float64)))])}

    def edge_end(self, a: dict[str, np.ndarray]) -> dict:
        """Forget a branch that was kept."""
        self.edges.pop(int(scalar(a, "scan")), None)
        return {}

    # ---- pairwise alignment ----------------------------------------------------------------------------------------------------

    def align_batch(self, a: dict[str, np.ndarray]) -> dict[str, np.ndarray]:
        """Fill the tables of affine-gap alignments (Gotoh: `src/seq/align.rs`) of one sequence, `centre`, against each of several, `bases`
        (all of them one after another, `lens` long each), in 32 bits, along the anti-diagonals of the tables for all the pairs at once.
        Every cell is made by the same additions, comparisons and choices, in the same order and with the same ties (M before X before Y), as
        the processor's, so the numbers are the same to the last bit; sequences that are shorter than the longest have columns on the right
        that nothing reads. The answer is what the program needs to read each alignment back as it does its own: for each pair the byte of
        every cell that says where it came from (`tm | tx << 2 | ty << 4`, for the three states), by row and then column of a table
        `longest + 1` wide, the three scores at the last column of each row (`n + 1` rows of three; the row 0 is not used), and at each
        column of the last row (`longest + 1` of three). Pairs need `n >= 1` and `lens >= 1`."""
        d = self.device
        centre = np.asarray(a["centre"], dtype=np.uint8)
        lens = np.asarray(a["lens"]).astype(np.int64)
        bases = np.asarray(a["bases"], dtype=np.uint8)
        free_ends = bool(int(scalar(a, "free_ends")))
        hit, mismatch, open_, extend = (np.float32(x) for x in a["scoring"])
        n, count, width = len(centre), len(lens), int(lens.max())
        if n < 1 or int(lens.min()) < 1:
            raise ValueError("a sequence with no base")

        def lower(x: np.ndarray) -> np.ndarray:
            return np.where((x >= 65) & (x <= 90), x + 32, x).astype(np.int16)

        padded = np.full((count, width), -1, dtype=np.int16)
        start = 0
        for k, m in enumerate(lens.tolist()):
            padded[k, :m] = lower(bases[start : start + m])
            start += m
        # the score of a gap run of each length, made with 32-bit operations in the order the processor makes them
        longest = max(n, width)
        runs = np.arange(longest + 1, dtype=np.float32)
        gap = np.zeros(longest + 1, dtype=np.float32) if free_ends else open_ + (runs - np.float32(1.0)) * extend
        gap = torch.from_numpy(gap).to(d)
        a_low = torch.from_numpy(lower(centre)).to(d)
        b_low = torch.from_numpy(padded).to(d)
        lens_t = torch.from_numpy(lens).to(d)
        pairs = torch.arange(count, device=d)
        rows = torch.arange(1, n + 1, device=d)
        neg = float("-inf")
        f32 = torch.float32
        match_t = torch.full((), float(hit), dtype=f32, device=d)
        mismatch_t = torch.full((), float(mismatch), dtype=f32, device=d)
        open_t = torch.full((), float(open_), dtype=f32, device=d)
        extend_t = torch.full((), float(extend), dtype=f32, device=d)

        def best3(m: torch.Tensor, x: torch.Tensor, y: torch.Tensor):
            # M against X, then the winner against Y: the same two choices as `best3` in `src/seq/align.rs`
            first = m >= x
            v = torch.where(first, m, x)
            t = torch.where(first, 0, 1)
            second = v >= y
            return torch.where(second, v, y), torch.where(second, t, 2)

        def table() -> list[torch.Tensor]:
            return [torch.full((count, n + 1), neg, dtype=f32, device=d) for _ in range(3)]

        trace = torch.zeros((count, (n + 1) * (width + 1)), dtype=torch.uint8, device=d)
        last_col = torch.full((count, n + 2, 3), neg, dtype=f32, device=d)
        last_row = torch.full((count, width + 1, 3), neg, dtype=f32, device=d)

        def record(diagonal: int, m: torch.Tensor, x: torch.Tensor, y: torch.Tensor) -> None:
            # the cell of each pair that stands at its last column on this diagonal (a spare place, `n + 1`, for the pairs with none)
            row = diagonal - lens_t
            here = (row >= 1) & (row <= n)
            at = row.clamp(0, n)[:, None]
            found = torch.stack([m.gather(1, at)[:, 0], x.gather(1, at)[:, 0], y.gather(1, at)[:, 0]], dim=1)
            last_col[pairs, torch.where(here, row, n + 1)] = found
            if diagonal >= n and diagonal - n <= width:
                last_row[:, diagonal - n, :] = torch.stack([m[:, n], x[:, n], y[:, n]], dim=1)

        # the diagonals 0 and 1: the corner, and the first cell of the first row and of the first column
        two = table()
        two[0][:, 0] = 0.0
        one_back = table()
        one_back[2][:, 0] = gap[1]
        one_back[1][:, 1] = gap[1]
        record(1, *one_back)
        for diagonal in range(2, n + width + 1):
            pm, px, py = two
            qm, qx, qy = one_back
            j_before = (diagonal - rows - 1).clamp(0, width - 1)
            same = a_low[None, :] == b_low.index_select(1, j_before)
            s = torch.where(same, match_t, mismatch_t)
            bm, tm = best3(pm[:, :n], px[:, :n], py[:, :n])
            bx, tx = best3(qm[:, :n] + open_t, qx[:, :n] + extend_t, qy[:, :n] + open_t)
            by, ty = best3(qm[:, 1:] + open_t, qx[:, 1:] + open_t, qy[:, 1:] + extend_t)
            new = table()
            new[0][:, 1:] = bm + s
            new[1][:, 1:] = bx
            new[2][:, 1:] = by
            # the first row (its place 0) and, where the diagonal reaches it, the first column
            new[2][:, 0] = gap[diagonal] if diagonal <= width else neg
            if diagonal <= n:
                new[0][:, diagonal] = neg
                new[1][:, diagonal] = gap[diagonal]
                new[2][:, diagonal] = neg
            low, high = max(1, diagonal - width), min(n, diagonal - 1)
            if low <= high:
                # (the three fields do not overlap, so sums are the `tm | tx << 2 | ty << 4` of the processor; MPS has no shifts of integers)
                code = (tm + tx * 4 + ty * 16)[:, low - 1 : high].to(torch.uint8)
                row = rows[low - 1 : high]
                trace[:, row * (width + 1) + (diagonal - row)] = code
            record(diagonal, *new)
            two, one_back = one_back, new
        # an error of the device shows here, not as the fault of the next request
        synchronize(self.device)
        return {
            "trace": trace.cpu().numpy().reshape(-1),
            "last_col": last_col[:, : n + 1, :].cpu().numpy().reshape(-1),
            "last_row": last_row.cpu().numpy().reshape(-1),
        }

    # ---- neighbour-joining -----------------------------------------------------------------------------------------------------

    def neighbor_joining(self, a: dict[str, np.ndarray]) -> dict[str, np.ndarray]:
        """The joins of neighbour-joining over the `n` by `n` distances `dist`, as `src/seq/tree.rs` makes them, in 64 bits: for each join the
        positions in the current list of clusters of the two joined (the new cluster goes to the end of the list) and their branch lengths, and
        the distances among the three clusters that are left. A joined cluster takes one of the two places of the matrix and the other is
        zeroed and ignored, so a join is a few passes over the matrix, not a copy of it; and when under 60 % of the places are still alive
        the matrix is cut down to them, so that the passes get shorter as the clusters get fewer. Where several pairs tie to the last digit,
        the pair earliest in the list is joined, as on the processor (whose sums are added in another order, so a tie to the last digit on
        one side can be a difference of one in the last on the other)."""
        if self.device.type == "mps":
            raise RuntimeError("this device has no 64-bit floats")
        n = int(scalar(a, "n"))
        d = self.device
        inf = float("inf")
        dist = torch.from_numpy(a["dist"].reshape(n, n).astype(np.float64)).to(d)
        buffer = torch.empty_like(dist)
        size = n
        active = torch.ones(size, dtype=torch.bool, device=d)
        # a place's rank is where its cluster stands in the list: lower is earlier
        rank = torch.arange(size, device=d)
        later = size
        order = list(range(size))  # places, earliest cluster first
        first, second, short, long = [], [], [], []
        for _ in range(n - 3):
            m = len(order)
            if 16 < m <= 0.6 * size:
                keep = torch.tensor(order, device=d)
                dist = dist.index_select(0, keep).index_select(1, keep).contiguous()
                buffer = torch.empty_like(dist)
                size = m
                active = torch.ones(size, dtype=torch.bool, device=d)
                rank = torch.arange(size, device=d)
                later = size
                order = list(range(size))
            r = dist.sum(1)
            # a place that is not alive has row sum minus infinity, which makes every q that involves it plus infinity
            dead = r.masked_fill(~active, -inf)
            torch.mul(dist, m - 2.0, out=buffer)
            buffer.sub_(dead[:, None]).sub_(dead[None, :])
            # only a pair whose first place is the earlier cluster counts, as the processor's q (m - 2) d - r_i - r_j is made for the earlier
            # one as i: the other order can differ in the last digit (the subtractions are not associative), and the diagonal is out too
            buffer.masked_fill_(rank[:, None] >= rank[None, :], inf)
            smallest = buffer.min()
            tied = buffer == smallest
            # of the pairs that tie, the one whose earlier cluster is earliest, and then whose later one is: the rows that have a tie, the
            # earliest of them, and in its row the earliest tie
            rows = tied.any(1)
            big = later + 1
            x = torch.argmin(torch.where(rows, rank, big))
            y = torch.argmin(torch.where(tied[x], rank, big))
            x, y = torch.stack([x, y]).tolist()
            px, py = order.index(x), order.index(y)
            dij, ri, rj = torch.stack([dist[x, y], r[x], r[y]]).tolist()
            m_less = m - 2.0
            li = max(dij / 2.0 + (ri - rj) / (2.0 * m_less), 0.0)
            lj = max(dij - li, 0.0)
            active[y] = False
            to_new = ((dist[x] + dist[y] - dij) / 2.0).masked_fill(~active, 0.0)
            to_new[x] = 0.0
            dist[x, :] = to_new
            dist[:, x] = to_new
            dist[y, :] = 0.0
            dist[:, y] = 0.0
            rank[x] = later
            later += 1
            order.remove(x)
            order.remove(y)
            order.append(x)
            first.append(px)
            second.append(py)
            short.append(li)
            long.append(lj)
        last = dist[order][:, order].cpu().numpy().astype(np.float64)
        return {
            "first": np.array(first, dtype="<u4"),
            "second": np.array(second, dtype="<u4"),
            "li": np.array(short, dtype=np.float64),
            "lj": np.array(long, dtype=np.float64),
            "last": last.reshape(-1),
        }

    def release(self, a: dict[str, np.ndarray]) -> dict:
        """Forget a set of columns. Its memory goes back to PyTorch's cache, to be used again; it goes back to the driver only when this
        process ends."""
        self.sets.pop(int(scalar(a, "id")), None)
        return {}


def main() -> None:
    device, kind, count = pick_device(sys.argv[1] if len(sys.argv) > 1 else "auto")
    if device.type == "cuda":
        torch.cuda.set_device(device)
    worker = Worker(device, kind)
    name, memory = describe(device, kind)
    write_message(
        "ok",
        {
            "count": np.array([count], dtype=np.uint32),
            "name": name,
            "kind": kind,
            "memory": np.array([float(memory)]),
            "torch": torch.__version__,
            "cuda": str(torch.version.cuda),
            "arches": " ".join(torch.cuda.get_arch_list()) if torch.cuda.is_available() else "",
        },
    )
    if count == 0:
        return
    handlers = {
        "speed": worker.speed,
        "load_distances": worker.load_distances,
        "counts": worker.counts,
        "load_likelihood": worker.load_likelihood,
        "likelihood": worker.likelihood,
        "edge_begin": worker.edge_begin,
        "edge_eval": worker.edge_eval,
        "edge_end": worker.edge_end,
        "align_batch": worker.align_batch,
        "neighbor_joining": worker.neighbor_joining,
        "release": worker.release,
    }
    while True:
        try:
            op, arrays = read_message()
        except EOFError:
            return
        reply, error, fatal = {}, None, False
        try:
            reply = handlers[op](arrays)
        except Exception as e:  # reported to the program, which goes on with the processor; this worker stays up for the other work
            traceback.print_exc()
            error = f"{type(e).__name__}: {e}"
            # out of memory can be tried again after the cache is emptied; another error of the card (a failed kernel, a lost device)
            # leaves it unusable, so this worker says so once and ends
            fatal = device.type == "cuda" and not isinstance(e, torch.cuda.OutOfMemoryError) and "CUDA" in str(e)
        # out of the `except`, so that the tensors of the request that failed are free by now
        if error is None:
            write_message("ok", reply)
            continue
        write_message("error", {"text": error})
        if device.type == "cuda":
            try:
                torch.cuda.empty_cache()
            except Exception:
                fatal = True
        if fatal:
            return


if __name__ == "__main__":
    main()
