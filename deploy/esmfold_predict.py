#!/usr/bin/env python3
"""Fold one protein with ESMFold on the GPU that the HelixLoom server gave this job.

The server runs a folding program as `PROGRAM input.fasta output.pdb` (see
`--fold-command` in DEPLOY.md). This is one such program. The server has already
limited the job to a single GPU through CUDA_VISIBLE_DEVICES, so the GPU seen here as
"cuda" is that one and no other.

Set-up on the workstation (the model weights are several gigabytes and are downloaded
the first time this runs):

    python3 -m venv /opt/helixloom/esmfold-env
    /opt/helixloom/esmfold-env/bin/pip install torch transformers accelerate

Written from the transformers documentation and source (`EsmForProteinFolding.infer_pdb`
takes one sequence and gives back the text of a PDB file, with the model's confidence
in the temperature column) but NOT run where it was written: there was no NVIDIA GPU
there. Try it by hand first:

    printf '>q\\nMQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG\\n' > q.fasta
    /opt/helixloom/esmfold-env/bin/python esmfold_predict.py q.fasta q.pdb
"""

import sys

import torch
from transformers import EsmForProteinFolding

MODEL = "facebook/esmfold_v1"
# Smaller chunks use less GPU memory and take longer: lower this if a long protein
# runs out of memory.
CHUNK_SIZE = 64


def read_sequence(path: str) -> str:
    """The residues of the first record of a FASTA file, as capital letters."""
    letters: list[str] = []
    with open(path, encoding="utf-8") as fasta:
        for line in fasta:
            if line.startswith(">"):
                if letters:
                    break
                continue
            letters.append(line.strip().upper())
    sequence = "".join(letters)
    if not sequence:
        sys.exit(f"{path} holds no sequence.")
    return sequence


def main() -> None:
    if len(sys.argv) != 3:
        sys.exit("usage: esmfold_predict.py input.fasta output.pdb")
    sequence = read_sequence(sys.argv[1])
    if not torch.cuda.is_available():
        sys.exit("No GPU is visible to this program, and folding on the processor is too slow.")
    model = EsmForProteinFolding.from_pretrained(MODEL, low_cpu_mem_usage=True)
    model = model.cuda()
    model.esm = model.esm.half()
    torch.backends.cuda.matmul.allow_tf32 = True
    model.trunk.set_chunk_size(CHUNK_SIZE)
    with torch.no_grad():
        pdb = model.infer_pdb(sequence)
    with open(sys.argv[2], "w", encoding="utf-8") as out:
        out.write(pdb)


if __name__ == "__main__":
    main()
