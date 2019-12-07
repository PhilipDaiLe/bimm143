Class11
================
Philip Dai Le
11/5/2019

\#\#PDB database practicce for biomolecular structure data \>**Q1**:
Download a CSV file from the PDB site (accessible from “Analyze” -\>
“PDB Statistics” \> “by Experimental Method and Molecular Type”.

> Move this CSV file into your RStudio project and determine the
> percentage of structures solved by X-Ray and Electron Microscopy.

Downloaded CSV file as “Data Export Summary.csv”

``` r
#reading in csv file
data<- read.csv("Data Export Summary.csv")
data
```

    ##   Experimental.Method Proteins Nucleic.Acids Protein.NA.Complex Other  Total
    ## 1               X-Ray   131278          2059               6759     8 140104
    ## 2                 NMR    11235          1303                261     8  12807
    ## 3 Electron Microscopy     2899            32                999     0   3930
    ## 4               Other      280             4                  6    13    303
    ## 5        Multi Method      144             5                  2     1    152

> Also can you determine what proportion of structures are protein? Aim
> to have a rendered GitHub document with working code that yields your
> answers.

``` r
sum(data$Total) #displays total of all
```

    ## [1] 157296

``` r
sum(data$Proteins) # displays total of protein
```

    ## [1] 145836

``` r
(sum(data$Proteins)/sum(data$Total))*100#displays proportion that is protein
```

    ## [1] 92.71437

``` r
(data$Total/sum(data$Total))*100 #proportion of entries for each method
```

    ## [1] 89.0702879  8.1419744  2.4984742  0.1926305  0.0966331

``` r
#two significant figs
round((sum(data$Proteins)/sum(data$Total))*100,2)
```

    ## [1] 92.71

> **Q2**: Type HIV in the PDB website search box on the home page and
> determine how many HIV-1 protease structures are in the current PDB?

\#\#HIV-1 PDB \>1.2 The PDB format Now download the “PDB File” for the
HIV-1 protease structure with the PDB identifier **1HSG**. On the
website you can “Display” the contents of this “PDB format” file.
Alternatively, you can examine the contents of your downloaded file in a
suitable text editor.

> Q3: Water molecules normally have 3 atoms. Why do we see just one atom
> per water molecule in this structure?

Due to the experimental method used to capture the image, hydrogen atoms
are too small to be depicted in the resolution.

> Q4: There is a conserved water molecule in the binding site. Can you
> identify this water molecule? What residue number does this water
> molecule have (see note below)?

Here we are reading in “1hsg.pdb” to select the protein component and
output a protein-only file. We’re then doing the same for the drug
ligand

``` r
library(bio3d)
pdb<-read.pdb("1hsg.pdb")
ligandSelect<- atom.select(pdb, "ligand", value=TRUE)
write.pdb(pdb)
#trim.pdb(pdb)
pdb #HOH is water, MK1 is ligand
```

    ## 
    ##  Call:  read.pdb(file = "1hsg.pdb")
    ## 
    ##    Total Models#: 1
    ##      Total Atoms#: 1686,  XYZs#: 5058  Chains#: 2  (values: A B)
    ## 
    ##      Protein Atoms#: 1514  (residues/Calpha atoms#: 198)
    ##      Nucleic acid Atoms#: 0  (residues/phosphate atoms#: 0)
    ## 
    ##      Non-protein/nucleic Atoms#: 172  (residues: 128)
    ##      Non-protein/nucleic resid values: [ HOH (127), MK1 (1) ]
    ## 
    ##    Protein sequence:
    ##       PQITLWQRPLVTIKIGGQLKEALLDTGADDTVLEEMSLPGRWKPKMIGGIGGFIKVRQYD
    ##       QILIEICGHKAIGTVLVGPTPVNIIGRNLLTQIGCTLNFPQITLWQRPLVTIKIGGQLKE
    ##       ALLDTGADDTVLEEMSLPGRWKPKMIGGIGGFIKVRQYDQILIEICGHKAIGTVLVGPTP
    ##       VNIIGRNLLTQIGCTLNF
    ## 
    ## + attr: atom, xyz, seqres, helix, sheet,
    ##         calpha, remark, call

``` r
?atom.select.pdb
```

    ## starting httpd help server ... done

``` r
write.pdb(ligandSelect, file="HSG1ligandSelect.pdb")
proteinSelect<-atom.select(pdb, "protein", value=TRUE)
write.pdb(proteinSelect, file="HSG1proteinSelect.pdb ")
```
