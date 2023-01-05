#!/bin/bash
./create_promoters_bed.py --gtf /data/CCBR_Pipeliner/db/PipeDB/Indices/GTFs/hg19/gencode.v42lift37.basic.annotation.gtf --out hg19_42.promoters.bed
./create_promoters_bed.py --gtf /data/CCBR_Pipeliner/db/PipeDB/Indices/GTFs/hg38/gencode.v42.primary_assembly.annotation.gtf --out hg38_42.promoters.bed
./create_promoters_bed.py --gtf /data/CCBR_Pipeliner/db/PipeDB/Indices/GTFs/mm10/gencode.vM25.annotation.gtf --out mm10_M25.promoters.bed
./create_promoters_bed.py --gtf /data/CCBR_Pipeliner/db/PipeDB/Indices/mmul10/genes.gtf --out mmul10.promoters.bed
./create_promoters_bed.py --gtf /data/CCBR_Pipeliner/db/PipeDB/Indices/bosTau9/genes.gtf --out bosTau9.promoters.bed
