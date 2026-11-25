# 🧬 Guppy Basecalling

## 📂 Input Directory Structure

Each input directory should follow this structure:

```
sample1/
├── fast5_pass/
├── fast5_fail/
├── fast5_skip/
```

---

## 🚀 Usage

### 1. Basecalling

```bash
python guppy_basecall.py basecall   --input_dirs /path/to/sample1 /path/to/sample2   --flowcell FLO-MIN106   --kit SQK-RNA002   --callers 4   --threads 8
```

#### Parameters

| Argument | Description | Default |
|----------|-------------|--------|
| `--input_dirs` | One or more input directories | Required |
| `--flowcell` | Flowcell type | FLO-MIN106 |
| `--kit` | Sequencing kit | SQK-RNA002 |
| `--callers` | Number of parallel basecallers | 2 |
| `--threads` | CPU threads per caller | 5 |

---

### 2. Merge FASTQ Files

```bash
python guppy_basecall.py merge   --input_dirs /path/to/sample1 /path/to/sample2
```

---

## 📤 Output Structure

After running the pipeline:

```
sample1/
├── guppy_out/
│   ├── pass/
│   ├── fail/
│   ├── skip/
│   └── (logs)
├── merge_guppy_outputs/
│   ├── All_pass.fastq
│   ├── All_fail.fastq
│   └── All_merged.fastq
```

---

## 🔁 Workflow Overview

1. Basecalling with Guppy (`fast5 → fastq`)
2. Merge FASTQ files:
   - Merge pass reads
   - Merge fail reads
   - Combine into one final FASTQ file

---

## 📌 Notes

- The script assumes Guppy is available in your system `PATH`
- Logs from basecalling are saved in each output folder (`guppy.log`)

---

