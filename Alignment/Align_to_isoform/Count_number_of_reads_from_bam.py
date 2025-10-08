import pysam
from collections import defaultdict

# 📥 Input and 📤 Output files
bam_file = "/data/groups/anko/mohamed/oxford_nanopore/Results/Nanopore/3_SRSF3_DIRCLIP/Alignment/bambu/mm9_transcriptome/CLIP_cDNA_rep1/gencode/filt_q8.bam"
gtf_file = "/data/groups/anko/mohamed/oxford_nanopore/genome/mm9/Annotation/gencode.vM1.annotation.gtf"
output_file = "/data/groups/anko/mohamed/oxford_nanopore/Results/Nanopore/3_SRSF3_DIRCLIP/Alignment/bambu/mm9_transcriptome/CLIP_cDNA_rep1/gencode/filt_q8_read_counts_per_transcript.tsv"

# -----------------------------------------------------------------------------
# 🧬 Step 1: Parse GTF file to create a transcript → gene name map
# -----------------------------------------------------------------------------
gene_id = "gene_name"
transcript_to_gene = {}

with open(gtf_file, 'r') as gtf:
    for line in gtf:
        if line.startswith("#"):
            continue
        fields = line.strip().split('\t')
        if fields[2] != "transcript":
            continue

        attributes = fields[8]
        # Extract key-value pairs from attribute string
        attrs = {
            kv.strip().split(' ')[0]: kv.strip().split(' ')[1].strip('"')
            for kv in attributes.strip(';').split('; ') if ' ' in kv
        }

        transcript_id = attrs.get("transcript_id")
        gene_name = attrs.get(gene_id, "NA")

        if transcript_id:
            transcript_to_gene[transcript_id] = gene_name

# -----------------------------------------------------------------------------
# 📊 Step 2: Count mapped reads per transcript
# -----------------------------------------------------------------------------
bam = pysam.AlignmentFile(bam_file, "rb")
transcript_counts = defaultdict(int)

for read in bam.fetch(until_eof=True):
    if not read.is_unmapped:
        transcript = read.reference_name

        # Strip "Chr_" prefix if present
        if transcript.startswith("Chr_"):
            transcript = transcript.replace("Chr_", "", 1)

        transcript_counts[transcript] += 1

bam.close()

# -----------------------------------------------------------------------------
# 📄 Step 3: Write transcript read counts and gene names to TSV output
# -----------------------------------------------------------------------------
with open(output_file, "w") as f:
    f.write("Isoform_id\tGene_id\tRead_counts\n")
    for transcript, count in sorted(transcript_counts.items()):
        gene = transcript_to_gene.get(transcript, "NA")
        f.write(f"{transcript}\t{gene}\t{count}\n")

print(f"✅ Saved read counts with gene names to: {output_file}")
