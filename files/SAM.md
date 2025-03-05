A **SAM (Sequence Alignment/Map) file** is a text-based format used to store alignments of sequencing reads to a reference genome. It is widely used in bioinformatics, particularly in next-generation sequencing (NGS) workflows.

---

### **Structure of a SAM File**
A SAM file consists of two main sections:

1. **Header Section (Optional)**
   - Begins with `@`
   - Contains metadata about the alignment, reference genome, and sequencing run.
   - Example:
     ```
     @HD   VN:1.6  SO:coordinate
     @SQ   SN:chr1 LN:248956422
     ```

2. **Alignment Section (Mandatory)**
   - Contains tab-separated fields describing each read’s alignment.
   - Example:
     ```
     read1   99   chr1   100   60   50M   =   200   100   AGCT...   IIIII...
     ```

---

### **Fields in the Alignment Section**
Each alignment entry contains **11 mandatory fields**:

| Column | Field Name        | Description |
|--------|------------------|-------------|
| 1      | QNAME (Query Name) | Read identifier (from FASTQ) |
| 2      | FLAG (Bitwise Flag) | Alignment details (e.g., paired-end, properly mapped) |
| 3      | RNAME (Reference Name) | Chromosome or contig name |
| 4      | POS (Position) | 1-based leftmost mapping position |
| 5      | MAPQ (Mapping Quality) | Confidence score of alignment |
| 6      | CIGAR (CIGAR String) | Representation of match, insertions, deletions |
| 7      | RNEXT (Reference Name of Mate) | Chromosome of mate read (for paired-end) |
| 8      | PNEXT (Position of Mate) | Position of mate read |
| 9      | TLEN (Template Length) | Insert size for paired-end reads |
| 10     | SEQ (Sequence) | Read sequence (from FASTQ) |
| 11     | QUAL (Quality Scores) | ASCII-encoded Phred quality |

---

### **Example SAM Entry**
```
SRR1234567  99  chr1  100  60  50M  =  200  100  AGCTTAGCTA...  IIIIIIIII...
```
- **`99`** → The FLAG value indicating this is a paired-end read.
- **`chr1`** → Read is mapped to chromosome 1.
- **`100`** → Read starts at position 100.
- **`60`** → Mapping quality score.
- **`50M`** → CIGAR string indicating a perfect 50-base match.
- **`AGCTTAGCTA...`** → Actual read sequence.
- **`IIIIIIII...`** → Quality scores.

---

### **Understanding the FLAG Field**
The **FLAG** field is a bitwise representation of alignment properties. Some common values:
- `0` → Unpaired read, mapped
- `16` → Read is mapped on the reverse strand
- `99` → Paired read, properly mapped
- `147` → Mate of a paired read, reverse strand

To decode FLAG values, tools like `samtools` (`samtools flagstat`) can be used.

---

### **Working with SAM Files**
- **Convert SAM to BAM (Binary Alignment/Map):**  
  ```
  samtools view -bS input.sam > output.bam
  ```
- **Sort BAM file:**  
  ```
  samtools sort output.bam -o output.sorted.bam
  ```
- **Index BAM file:**  
  ```
  samtools index output.sorted.bam
  ```
- **View alignments:**  
  ```
  samtools tview output.sorted.bam reference.fasta
  ```

---
