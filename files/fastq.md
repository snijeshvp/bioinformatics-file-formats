A **FASTQ file** is a widely used format in bioinformatics for storing raw sequencing data from high-throughput sequencing platforms like Illumina, PacBio, and Oxford Nanopore. It contains both nucleotide sequences and their associated quality scores.

### **FASTQ File Structure**
Each sequence in a FASTQ file consists of **four lines**:
1. **Header Line (Identifier)**
   - Starts with `@`
   - Contains information about the sequencing run, lane, and read number.
   - Example:  
     ```
     @SEQ_ID
     ```

2. **Sequence Line**
   - Contains the nucleotide sequence (A, T, C, G, or N for unknown bases).
   - Example:  
     ```
     GATTTGGGGTTTCCCAGTCACGAC
     ```

3. **Quality Header Line**
   - Starts with `+`
   - May optionally repeat the sequence identifier.


4. **Quality Score Line**
   - Encodes per-base quality scores using ASCII characters.
   - Higher ASCII values indicate better base call confidence.
   - Example:  
     ```
     !''*((((***+))%%%++)(%%%%).1
     ```

### **Quality Score Encoding**
The quality score is represented using **Phred scores**, which are encoded as ASCII characters. The Phred score is calculated as:

$$
Q = -10 \times \log_{10}(P)
$$

where **P** is the probability of an incorrect base call.

- Common Phred scoring schemes:
  - **Sanger (Phred+33)**: ASCII 33–73 (Used in most modern datasets)
  - **Illumina 1.3+ (Phred+64)**: ASCII 64–104 (Older format)

### **Example FASTQ Entry**
```
@SEQ_ID
GATTTGGGGTTTCCCAGTCACGAC
+
!''*((((***+))%%%++)(%%%%).1
```
