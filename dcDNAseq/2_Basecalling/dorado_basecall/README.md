# 🧬 Dorado Basecalling
---

## 📦 Requirements

- Python ≥ 3.7
- Oxford Nanopore **Dorado basecaller**
- samtools
- pigz
- Linux/Unix environment

---

### Run full workflow
```bash
python dorado_basecaller.py run   --input_dirs /path/to/pod5_dir   --tmp_dir /path/to/tmp   --output_dir /path/to/output   --model /path/to/dorado_model   --file_type .pod5   --batch_size 4
```

---

## 📤 Output Structure

After running the pipeline:

```
output/
├── *.bam
├── merged/
│   ├── All_merged.bam
│   ├── All_merged.fastq.gz
│   └── summary.tsv
```

---

## 🔁 Workflow Overview

1. Batch basecalling (pod5/fast5 → BAM)
2. Merge BAM files
3. Generate summary statistics
4. Convert BAM → FASTQ (gzipped)

---

## 📌 Notes

- Ensure `dorado`, and `samtools` are in your PATH

---

