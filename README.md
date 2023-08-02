# UDP_RNAseq

## Purpose

This document should aid the user in accessing Argon HPC, retreiving RNA-seq data from IIHG, and generating counts for differential expression analysis

## Accessing Argon

Argon access must be applied for through ITS. This can be done through [this link](https://research.its.uiowa.edu/our-services/computing-services/argon-high-performance-computing-hpc).

Once access is granted, Argon can be accessed through SecureCRT, available for download [here](https://helpdesk.its.uiowa.edu/software/download/securecrt/). Be sure to download both the program file and the license file. The license information will be requested the first time you run the program. Simply copy the contents of the license file and paste them into the indicated box.

### SecureCRT

In SecureCRT, to access Argon, you must set up the ssh connection. To do this, go to File > Connect.Under Session Manager, click the + sign to open the New Session Wizard.

1. Set the Protocol to SSH2, click next.
2. Hostname: argon.hpc.uiowa.edu
3. Port can be left as the default (22). If you are working from an off-campus location, set this to 40.
4. Firewall: none
5. Username: Use your HawkID.
6. Click Next, then Finish.

The system will request that you enter your password - it's the same as your UIowa or Healthcare login. Two factor authentication must be set up through DUO.

## The Linux Command Line

To use Argon, you must be familiar with the linux command line. Some useful commands are:

* `cd {location}` Allows you to change your directory. This allows you to move from folder to folder.
    + `cd ..` Moves you up to the parent directory.
    + `cd ~` Moves you to your home directory. This can also be acheived using `cd /Users/HAWKID`
    + `cd /Shared/lss_chndr` Moves you to the LSS drive called lss_chndr.
* `ls` Allows you to see what directories and files are available in your current directory.
* `cp {file} {destination}` Allows you to copy a file to a new location.
* `mv {file} {destination}` Allows you to move a file to a new location.
* `pwd` Shows you your current location.
* `mkdir {name}` Create a new directory (folder).
* `rm {file}` Delete a file.
    + `rm -r {directory}` To delete a folder, you must specify `-r` for recursive removal.

Note that LSS drives are accessible outside of Argon by mapping them to your computer. This means you can manipulate the contents of your LSS drive both through the Argon command line interface and the Windows/Apple GUI. Any files moved from your home directory `/Users/HAWKID` to the LSS `/Shared/lss_chndr` become accessible to your computer.

## Downloading Data from IIHG

You will receive an email from IIHG with a link to your data. In Argon, navigate to a directory in which you would like to work.

1. `mkdir Downloads`
    i) This creates a new directory (folder) to put our downloaded data in.
2. `cd Downloads`
    i) Move into that new directory.
3. `wget -r -np {link}`
    i) Paste your link from IIHG in place of {link}. This will begin to download all of the sample data from the link.
4. `cat Lane1/SampleName_R1.fastq.gz Lane2/SampleName_R1.fastq.gz > SampleName_R1.fastq.gz`
    i) Concatenate the fastq files from lane 1 and lane 2. Match based on sample name and replicate number. Save those as a new file that has the sample name and replicate number.

## Create an index

This will use Salmon. The guide for Salmon can be found [here](https://combine-lab.github.io/salmon/getting_started/). Files will be generated from this, so you may want to move to a useful directory or make a new one and move into it.

1. `module load stack/2022.2`
    i) This is necessary for the next step to work correctly.
2. `module load salmon`
    i) Salmon is already available on Argon; we simply need to load it to make it available for use.
3. `curl https://ftp.ensembl.org/pub/release-109/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz -o GRCh38_cdna.fa.gz`
    i) This downloads a cDNA reference from Ensembl. The link may need to be updated in the future.
4. `salmon index -t GRCh38_cdna.fa.gz -i GRCh38_index`
    i) This creates an index called GRCh38_index to map RNA-seq reads to for quantification.

## Quantify Samples

We will use our index we created as well as our concatenated fastq files for this step. It's important to include the location for each when referencing them. I'll pretend they're in directories called `/Downloads/` and `/Index/` for this.

You will have multiple samples and running this takes a while. The easiest way to do this is to submit a job through Argon that can run without your oversight. To do this, first we need to make a job file. This can just be done in a text editor such as Notepad and saved to the LSS for you to access later through Argon. You can make a shell script similar to the one below:

```
#!/bin/bash

module load stack/2022.2
module load salmon

for f in SampleName1 SampleName2 SampleName3 SampleName4
do
salmon quant -i ~/Index/GRCh38_index -l A \
       	-1 ~/Downloads/${f}_R1.fastq.gz \
       	-2 ~/Downloads/${f}_R2.fastq.gz \
       	-p 8 --validateMappings -o quants/${f}_quant
done
```

This assumes you have four samples with two files each. For example you have SampleName1_R1.fastq.gz and SampleName1_R2.fastq.gz. You can add as many samples as you'd like, just be sure to update the names. Save this file as a .job file in the LSS.

Once you have that saved:
1. In your home directory in Argon (`cd ~`), create a new directory called quants (`mkdir quants`)
2. Navigate to the LSS (`cd /Shared/lss_chndr`).
3. `sed -i 's/\r$//' myscript.job`
4. `qsub myscript.job`
    i) This will run the script you created. Be sure to edit the command based on what you named your .job file.
5. `cd ~`
6. `ls`
    i) If the job has started (it may take a minute), you should see a file similar to myjob.job.e###
7. `cat myjob.job.e###`
    i) This is where the system messages related to the job will appear. If there are error codes, they will appear here and any log that Salmon creates will appear here as well.
    
Now it is a waiting game - Salmon seems to take at least 30 minutes per sample. Since we submitted this as a job, Argon can be closed and returned to later.

Congrats! You've completed transcript quantification!
