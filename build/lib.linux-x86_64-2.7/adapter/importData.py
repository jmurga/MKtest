import pyfaidx as px
import numpy as np
import pandas as pd

class MultiFasta:

	def __init__(self,file,codonTable,split,bins):
		self.file       = file
		self.codonTable = codonTable
		self.split      = split
		self.bins       = bins

	def degenerancy(self,data):
		
		#DEGENERANCY DICTIONARIES
		standardDict = {
		'TTT': '002', 'TTC': '002', 'TTA': '202', 'TTG': '202', 'TCT': '004', 'TCC': '004', 'TCA': '004', 'TCG': '004', 'TAT': '002', 'TAC': '002', 'TAA': '022', 'TAG': '002', 'TGT': '002', 'TGC': '002', 'TGA': '020', 'TGG': '000', 'CTT': '004', 'CTC': '004', 'CTA': '204', 'CTG': '204', 'CCT': '004', 'CCC': '004', 'CCA': '004', 'CCG': '004', 'CAT': '002', 'CAC': '002', 'CAA': '002', 'CAG': '002', 'CGT': '004', 'CGC': '004', 'CGA': '204', 'CGG': '204', 'ATT': '003', 'ATC': '003', 'ATA': '003', 'ATG': '000', 'ACT': '004', 'ACC': '004', 'ACA': '004', 'ACG': '004', 'AAT': '002', 'AAC': '002', 'AAA': '002', 'AAG': '002', 'AGT': '002', 'AGC': '002', 'AGA': '202', 'AGG': '202', 'GTT': '004', 'GTC': '004', 'GTA': '004', 'GTG': '004', 'GCT': '004', 'GCC': '004', 'GCA': '004', 'GCG': '004', 'GAT': '002', 'GAC': '002', 'GAA': '002', 'GAG': '002', 'GGT': '004', 'GGC': '004', 'GGA': '004', 'GGG': '004'}
		VertebrateMitochondrial = {

			'TTT': '002', 'TTC': '002', 'TTA': '202', 'TTG': '202', 'TCT': '004', 'TCC': '004', 'TCA': '004', 'TCG': '004', 'TAT': '002', 'TAC': '002', 'TAA': '002', 'TAG': '002', 'TGT': '002', 'TGC': '002', 'TGA': '002', 'TGG': '002', 'CTT': '004', 'CTC': '004', 'CTA': '204', 'CTG': '204', 'CCT': '004', 'CCC': '004', 'CCA': '004', 'CCG': '004', 'CAT': '002', 'CAC': '002', 'CAA': '002', 'CAG': '002', 'CGT': '004', 'CGC': '004', 'CGA': '004', 'CGG': '004', 'ATT': '002', 'ATC': '002', 'ATA': '002', 'ATG': '002', 'ACT': '004', 'ACC': '004', 'ACA': '004', 'ACG': '004', 'AAT': '002', 'AAC': '002', 'AAA': '002', 'AAG': '002', 'AGT': '002', 'AGC': '002', 'AGA': '002', 'AGG': '002', 'GTT': '004', 'GTC': '004', 'GTA': '004', 'GTG': '004', 'GCT': '004', 'GCC': '004', 'GCA': '004', 'GCG': '004', 'GAT': '002', 'GAC': '002', 'GAA': '002', 'GAG': '002', 'GGT': '004', 'GGC': '004', 'GGA': '004', 'GGG': '004'}
		YeastMitochondrial = {

			'TTT': '002', 'TTC': '002', 'TTA': '002', 'TTG': '002', 'TCT': '004', 'TCC': '004', 'TCA': '004', 'TCG': '004', 'TAT': '002', 'TAC': '002', 'TAA': '002', 'TAG': '002', 'TGT': '002', 'TGC': '002', 'TGA': '002', 'TGG': '002', 'CTT': '004', 'CTC': '004', 'CTA': '004', 'CTG': '004', 'CCT': '004', 'CCC': '004', 'CCA': '004', 'CCG': '004', 'CAT': '002', 'CAC': '002', 'CAA': '002', 'CAG': '002', 'CGT': '004', 'CGC': '004', 'CGA': '204', 'CGG': '204', 'ATT': '002', 'ATC': '002', 'ATA': '002', 'ATG': '002', 'ACT': '004', 'ACC': '004', 'ACA': '004', 'ACG': '004', 'AAT': '002', 'AAC': '002', 'AAA': '002', 'AAG': '002', 'AGT': '002', 'AGC': '002', 'AGA': '202', 'AGG': '202', 'GTT': '004', 'GTC': '004', 'GTA': '004', 'GTG': '004', 'GCT': '004', 'GCC': '004', 'GCA': '004', 'GCG': '004', 'GAT': '002', 'GAC': '002', 'GAA': '002', 'GAG': '002', 'GGT': '004', 'GGC': '004', 'GGA': '004', 'GGG': '004'}
		MoldMitochondrialProtozoanMitochondrialCoelenterateMitochondrialMycoplasmaSpiroplasma = {

			'TTT': '002', 'TTC': '002', 'TTA': '202', 'TTG': '202', 'TCT': '004', 'TCC': '004', 'TCA': '004', 'TCG': '004', 'TAT': '002', 'TAC': '002', 'TAA': '002', 'TAG': '002', 'TGT': '002', 'TGC': '002', 'TGA': '002', 'TGG': '002', 'CTT': '004', 'CTC': '004', 'CTA': '204', 'CTG': '204', 'CCT': '004', 'CCC': '004', 'CCA': '004', 'CCG': '004', 'CAT': '002', 'CAC': '002', 'CAA': '002', 'CAG': '002', 'CGT': '004', 'CGC': '004', 'CGA': '204', 'CGG': '204', 'ATT': '003', 'ATC': '003', 'ATA': '003', 'ATG': '000', 'ACT': '004', 'ACC': '004', 'ACA': '004', 'ACG': '004', 'AAT': '002', 'AAC': '002', 'AAA': '002', 'AAG': '002', 'AGT': '002', 'AGC': '002', 'AGA': '202', 'AGG': '202', 'GTT': '004', 'GTC': '004', 'GTA': '004', 'GTG': '004', 'GCT': '004', 'GCC': '004', 'GCA': '004', 'GCG': '004', 'GAT': '002', 'GAC': '002', 'GAA': '002', 'GAG': '002', 'GGT': '004', 'GGC': '004', 'GGA': '004', 'GGG': '004'}
		InvertebrateMitochondrial = {

			'TTT': '002', 'TTC': '002', 'TTA': '202', 'TTG': '202', 'TCT': '004', 'TCC': '004', 'TCA': '004', 'TCG': '004', 'TAT': '002', 'TAC': '002', 'TAA': '002', 'TAG': '002', 'TGT': '002', 'TGC': '002', 'TGA': '002', 'TGG': '002', 'CTT': '004', 'CTC': '004', 'CTA': '204', 'CTG': '204', 'CCT': '004', 'CCC': '004', 'CCA': '004', 'CCG': '004', 'CAT': '002', 'CAC': '002', 'CAA': '002', 'CAG': '002', 'CGT': '004', 'CGC': '004', 'CGA': '004', 'CGG': '004', 'ATT': '002', 'ATC': '002', 'ATA': '002', 'ATG': '002', 'ACT': '004', 'ACC': '004', 'ACA': '004', 'ACG': '004', 'AAT': '002', 'AAC': '002', 'AAA': '002', 'AAG': '002', 'AGT': '004', 'AGC': '004', 'AGA': '004', 'AGG': '004', 'GTT': '004', 'GTC': '004', 'GTA': '004', 'GTG': '004', 'GCT': '004', 'GCC': '004', 'GCA': '004', 'GCG': '004', 'GAT': '002', 'GAC': '002', 'GAA': '002', 'GAG': '002', 'GGT': '004', 'GGC': '004', 'GGA': '004', 'GGG': '004'}
		CiliateNuclearDasycladaceanNuclearHexamitaNuclear = {

			'TTT': '002', 'TTC': '002', 'TTA': '202', 'TTG': '202', 'TCT': '004', 'TCC': '004', 'TCA': '004', 'TCG': '004', 'TAT': '002', 'TAC': '002', 'TAA': '202', 'TAG': '202', 'TGT': '002', 'TGC': '002', 'TGA': '000', 'TGG': '000', 'CTT': '004', 'CTC': '004', 'CTA': '204', 'CTG': '204', 'CCT': '004', 'CCC': '004', 'CCA': '004', 'CCG': '004', 'CAT': '002', 'CAC': '002', 'CAA': '202', 'CAG': '202', 'CGT': '004', 'CGC': '004', 'CGA': '204', 'CGG': '204', 'ATT': '003', 'ATC': '003', 'ATA': '003', 'ATG': '000', 'ACT': '004', 'ACC': '004', 'ACA': '004', 'ACG': '004', 'AAT': '002', 'AAC': '002', 'AAA': '002', 'AAG': '002', 'AGT': '002', 'AGC': '002', 'AGA': '202', 'AGG': '202', 'GTT': '004', 'GTC': '004', 'GTA': '004', 'GTG': '004', 'GCT': '004', 'GCC': '004', 'GCA': '004', 'GCG': '004', 'GAT': '002', 'GAC': '002', 'GAA': '002', 'GAG': '002', 'GGT': '004', 'GGC': '004', 'GGA': '004', 'GGG': '004'}
		EchinodermMitochondrialFlatwormMitochondrial = {

			'TTT': '002', 'TTC': '002', 'TTA': '202', 'TTG': '202', 'TCT': '004', 'TCC': '004', 'TCA': '004', 'TCG': '004', 'TAT': '002', 'TAC': '002', 'TAA': '002', 'TAG': '002', 'TGT': '002', 'TGC': '002', 'TGA': '002', 'TGG': '002', 'CTT': '004', 'CTC': '004', 'CTA': '204', 'CTG': '204', 'CCT': '004', 'CCC': '004', 'CCA': '004', 'CCG': '004', 'CAT': '002', 'CAC': '002', 'CAA': '002', 'CAG': '002', 'CGT': '004', 'CGC': '004', 'CGA': '004', 'CGG': '004', 'ATT': '003', 'ATC': '003', 'ATA': '003', 'ATG': '000', 'ACT': '004', 'ACC': '004', 'ACA': '004', 'ACG': '004', 'AAT': '003', 'AAC': '003', 'AAA': '003', 'AAG': '000', 'AGT': '004', 'AGC': '004', 'AGA': '004', 'AGG': '004', 'GTT': '004', 'GTC': '004', 'GTA': '004', 'GTG': '004', 'GCT': '004', 'GCC': '004', 'GCA': '004', 'GCG': '004', 'GAT': '002', 'GAC': '002', 'GAA': '002', 'GAG': '002', 'GGT': '004', 'GGC': '004', 'GGA': '004', 'GGG': '004'}
		EuplotidNuclear = {

			'TTT': '002', 'TTC': '002', 'TTA': '202', 'TTG': '202', 'TCT': '004', 'TCC': '004', 'TCA': '004', 'TCG': '004', 'TAT': '002', 'TAC': '002', 'TAA': '002', 'TAG': '002', 'TGT': '003', 'TGC': '003', 'TGA': '003', 'TGG': '000', 'CTT': '004', 'CTC': '004', 'CTA': '204', 'CTG': '204', 'CCT': '004', 'CCC': '004', 'CCA': '004', 'CCG': '004', 'CAT': '002', 'CAC': '002', 'CAA': '002', 'CAG': '002', 'CGT': '004', 'CGC': '004', 'CGA': '204', 'CGG': '204', 'ATT': '003', 'ATC': '003', 'ATA': '003', 'ATG': '000', 'ACT': '004', 'ACC': '004', 'ACA': '004', 'ACG': '004', 'AAT': '002', 'AAC': '002', 'AAA': '002', 'AAG': '002', 'AGT': '002', 'AGC': '002', 'AGA': '202', 'AGG': '202', 'GTT': '004', 'GTC': '004', 'GTA': '004', 'GTG': '004', 'GCT': '004', 'GCC': '004', 'GCA': '004', 'GCG': '004', 'GAT': '002', 'GAC': '002', 'GAA': '002', 'GAG': '002', 'GGT': '004', 'GGC': '004', 'GGA': '004', 'GGG': '004'}
		BacterialArchaealandPlantPlastid = {

			'TTT': '002', 'TTC': '002', 'TTA': '202', 'TTG': '202', 'TCT': '004', 'TCC': '004', 'TCA': '004', 'TCG': '004', 'TAT': '002', 'TAC': '002', 'TAA': '022', 'TAG': '002', 'TGT': '002', 'TGC': '002', 'TGA': '020', 'TGG': '000', 'CTT': '004', 'CTC': '004', 'CTA': '204', 'CTG': '204', 'CCT': '004', 'CCC': '004', 'CCA': '004', 'CCG': '004', 'CAT': '002', 'CAC': '002', 'CAA': '002', 'CAG': '002', 'CGT': '004', 'CGC': '004', 'CGA': '204', 'CGG': '204', 'ATT': '003', 'ATC': '003', 'ATA': '003', 'ATG': '000', 'ACT': '004', 'ACC': '004', 'ACA': '004', 'ACG': '004', 'AAT': '002', 'AAC': '002', 'AAA': '002', 'AAG': '002', 'AGT': '002', 'AGC': '002', 'AGA': '202', 'AGG': '202', 'GTT': '004', 'GTC': '004', 'GTA': '004', 'GTG': '004', 'GCT': '004', 'GCC': '004', 'GCA': '004', 'GCG': '004', 'GAT': '002', 'GAC': '002', 'GAA': '002', 'GAG': '002', 'GGT': '004', 'GGC': '004', 'GGA': '004', 'GGG': '004'}
		AlternativeYeastNuclear = {

			'TTT': '002', 'TTC': '002', 'TTA': '202', 'TTG': '002', 'TCT': '004', 'TCC': '004', 'TCA': '004', 'TCG': '004', 'TAT': '002', 'TAC': '002', 'TAA': '022', 'TAG': '002', 'TGT': '002', 'TGC': '002', 'TGA': '020', 'TGG': '000', 'CTT': '003', 'CTC': '003', 'CTA': '203', 'CTG': '000', 'CCT': '004', 'CCC': '004', 'CCA': '004', 'CCG': '004', 'CAT': '002', 'CAC': '002', 'CAA': '002', 'CAG': '002', 'CGT': '004', 'CGC': '004', 'CGA': '204', 'CGG': '204', 'ATT': '003', 'ATC': '003', 'ATA': '003', 'ATG': '000', 'ACT': '004', 'ACC': '004', 'ACA': '004', 'ACG': '004', 'AAT': '002', 'AAC': '002', 'AAA': '002', 'AAG': '002', 'AGT': '002', 'AGC': '002', 'AGA': '202', 'AGG': '202', 'GTT': '004', 'GTC': '004', 'GTA': '004', 'GTG': '004', 'GCT': '004', 'GCC': '004', 'GCA': '004', 'GCG': '004', 'GAT': '002', 'GAC': '002', 'GAA': '002', 'GAG': '002', 'GGT': '004', 'GGC': '004', 'GGA': '004', 'GGG': '004'}
		AlternativeFlatwormMitochondrial = {

			'TTT': '002', 'TTC': '002', 'TTA': '202', 'TTG': '202', 'TCT': '004', 'TCC': '004', 'TCA': '004', 'TCG': '004', 'TAT': '003', 'TAC': '003', 'TAA': '003', 'TAG': '000', 'TGT': '002', 'TGC': '002', 'TGA': '002', 'TGG': '002', 'CTT': '004', 'CTC': '004', 'CTA': '204', 'CTG': '204', 'CCT': '004', 'CCC': '004', 'CCA': '004', 'CCG': '004', 'CAT': '002', 'CAC': '002', 'CAA': '002', 'CAG': '002', 'CGT': '004', 'CGC': '004', 'CGA': '004', 'CGG': '004', 'ATT': '003', 'ATC': '003', 'ATA': '003', 'ATG': '000', 'ACT': '004', 'ACC': '004', 'ACA': '004', 'ACG': '004', 'AAT': '003', 'AAC': '003', 'AAA': '003', 'AAG': '000', 'AGT': '004', 'AGC': '004', 'AGA': '004', 'AGG': '004', 'GTT': '004', 'GTC': '004', 'GTA': '004', 'GTG': '004', 'GCT': '004', 'GCC': '004', 'GCA': '004', 'GCG': '004', 'GAT': '002', 'GAC': '002', 'GAA': '002', 'GAG': '002', 'GGT': '004', 'GGC': '004', 'GGA': '004', 'GGG': '004'}
		BlepharismaMacronuclear = {

			'TTT': '002', 'TTC': '002', 'TTA': '202', 'TTG': '202', 'TCT': '004', 'TCC': '004', 'TCA': '004', 'TCG': '004', 'TAT': '002', 'TAC': '002', 'TAA': '020', 'TAG': '200', 'TGT': '002', 'TGC': '002', 'TGA': '020', 'TGG': '000', 'CTT': '004', 'CTC': '004', 'CTA': '204', 'CTG': '204', 'CCT': '004', 'CCC': '004', 'CCA': '004', 'CCG': '004', 'CAT': '002', 'CAC': '002', 'CAA': '002', 'CAG': '202', 'CGT': '004', 'CGC': '004', 'CGA': '204', 'CGG': '204', 'ATT': '003', 'ATC': '003', 'ATA': '003', 'ATG': '000', 'ACT': '004', 'ACC': '004', 'ACA': '004', 'ACG': '004', 'AAT': '002', 'AAC': '002', 'AAA': '002', 'AAG': '002', 'AGT': '002', 'AGC': '002', 'AGA': '202', 'AGG': '202', 'GTT': '004', 'GTC': '004', 'GTA': '004', 'GTG': '004', 'GCT': '004', 'GCC': '004', 'GCA': '004', 'GCG': '004', 'GAT': '002', 'GAC': '002', 'GAA': '002', 'GAG': '002', 'GGT': '004', 'GGC': '004', 'GGA': '004', 'GGG': '004'}
		ChlorophyceanMitochondrial = {

			'TTT': '002', 'TTC': '002', 'TTA': '202', 'TTG': '222', 'TCT': '004', 'TCC': '004', 'TCA': '004', 'TCG': '004', 'TAT': '002', 'TAC': '002', 'TAA': '020', 'TAG': '020', 'TGT': '002', 'TGC': '002', 'TGA': '020', 'TGG': '000', 'CTT': '004', 'CTC': '004', 'CTA': '204', 'CTG': '204', 'CCT': '004', 'CCC': '004', 'CCA': '004', 'CCG': '004', 'CAT': '002', 'CAC': '002', 'CAA': '002', 'CAG': '002', 'CGT': '004', 'CGC': '004', 'CGA': '204', 'CGG': '204', 'ATT': '003', 'ATC': '003', 'ATA': '003', 'ATG': '000', 'ACT': '004', 'ACC': '004', 'ACA': '004', 'ACG': '004', 'AAT': '002', 'AAC': '002', 'AAA': '002', 'AAG': '002', 'AGT': '002', 'AGC': '002', 'AGA': '202', 'AGG': '202', 'GTT': '004', 'GTC': '004', 'GTA': '004', 'GTG': '004', 'GCT': '004', 'GCC': '004', 'GCA': '004', 'GCG': '004', 'GAT': '002', 'GAC': '002', 'GAA': '002', 'GAG': '002', 'GGT': '004', 'GGC': '004', 'GGA': '004', 'GGG': '004'}
		TrematodeMitochondrial = {

			'TTT': '002', 'TTC': '002', 'TTA': '202', 'TTG': '202', 'TCT': '004', 'TCC': '004', 'TCA': '004', 'TCG': '004', 'TAT': '002', 'TAC': '002', 'TAA': '002', 'TAG': '002', 'TGT': '002', 'TGC': '002', 'TGA': '002', 'TGG': '002', 'CTT': '004', 'CTC': '004', 'CTA': '204', 'CTG': '204', 'CCT': '004', 'CCC': '004', 'CCA': '004', 'CCG': '004', 'CAT': '002', 'CAC': '002', 'CAA': '002', 'CAG': '002', 'CGT': '004', 'CGC': '004', 'CGA': '004', 'CGG': '004', 'ATT': '002', 'ATC': '002', 'ATA': '002', 'ATG': '002', 'ACT': '004', 'ACC': '004', 'ACA': '004', 'ACG': '004', 'AAT': '003', 'AAC': '003', 'AAA': '003', 'AAG': '000', 'AGT': '004', 'AGC': '004', 'AGA': '004', 'AGG': '004', 'GTT': '004', 'GTC': '004', 'GTA': '004', 'GTG': '004', 'GCT': '004', 'GCC': '004', 'GCA': '004', 'GCG': '004', 'GAT': '002', 'GAC': '002', 'GAA': '002', 'GAG': '002', 'GGT': '004', 'GGC': '004', 'GGA': '004', 'GGG': '004'}
		ScenedesmusobliquusMitochondrial = {

			'TTT': '002', 'TTC': '002', 'TTA': '202', 'TTG': '222', 'TCT': '003', 'TCC': '003', 'TCA': '030', 'TCG': '003', 'TAT': '002', 'TAC': '002', 'TAA': '030', 'TAG': '020', 'TGT': '002', 'TGC': '002', 'TGA': '030', 'TGG': '000', 'CTT': '004', 'CTC': '004', 'CTA': '204', 'CTG': '204', 'CCT': '004', 'CCC': '004', 'CCA': '004', 'CCG': '004', 'CAT': '002', 'CAC': '002', 'CAA': '002', 'CAG': '002', 'CGT': '004', 'CGC': '004', 'CGA': '204', 'CGG': '204', 'ATT': '003', 'ATC': '003', 'ATA': '003', 'ATG': '000', 'ACT': '004', 'ACC': '004', 'ACA': '004', 'ACG': '004', 'AAT': '002', 'AAC': '002', 'AAA': '002', 'AAG': '002', 'AGT': '002', 'AGC': '002', 'AGA': '202', 'AGG': '202', 'GTT': '004', 'GTC': '004', 'GTA': '004', 'GTG': '004', 'GCT': '004', 'GCC': '004', 'GCA': '004', 'GCG': '004', 'GAT': '002', 'GAC': '002', 'GAA': '002', 'GAG': '002', 'GGT': '004', 'GGC': '004', 'GGA': '004', 'GGG': '004'}
		ThraustochytriumMitochondrial = {

			'TTT': '002', 'TTC': '002', 'TTA': '030', 'TTG': '200', 'TCT': '004', 'TCC': '004', 'TCA': '004', 'TCG': '004', 'TAT': '002', 'TAC': '002', 'TAA': '032', 'TAG': '002', 'TGT': '002', 'TGC': '002', 'TGA': '030', 'TGG': '000', 'CTT': '004', 'CTC': '004', 'CTA': '004', 'CTG': '204', 'CCT': '004', 'CCC': '004', 'CCA': '004', 'CCG': '004', 'CAT': '002', 'CAC': '002', 'CAA': '002', 'CAG': '002', 'CGT': '004', 'CGC': '004', 'CGA': '204', 'CGG': '204', 'ATT': '003', 'ATC': '003', 'ATA': '003', 'ATG': '000', 'ACT': '004', 'ACC': '004', 'ACA': '004', 'ACG': '004', 'AAT': '002', 'AAC': '002', 'AAA': '002', 'AAG': '002', 'AGT': '002', 'AGC': '002', 'AGA': '202', 'AGG': '202', 'GTT': '004', 'GTC': '004', 'GTA': '004', 'GTG': '004', 'GCT': '004', 'GCC': '004', 'GCA': '004', 'GCG': '004', 'GAT': '002', 'GAC': '002', 'GAA': '002', 'GAG': '002', 'GGT': '004', 'GGC': '004', 'GGA': '004', 'GGG': '004'}
		PterobranchiaMitochondrial = {

			'TTT': '002', 'TTC': '002', 'TTA': '202', 'TTG': '202', 'TCT': '004', 'TCC': '004', 'TCA': '004', 'TCG': '004', 'TAT': '002', 'TAC': '002', 'TAA': '002', 'TAG': '002', 'TGT': '002', 'TGC': '002', 'TGA': '002', 'TGG': '002', 'CTT': '004', 'CTC': '004', 'CTA': '204', 'CTG': '204', 'CCT': '004', 'CCC': '004', 'CCA': '004', 'CCG': '004', 'CAT': '002', 'CAC': '002', 'CAA': '002', 'CAG': '002', 'CGT': '004', 'CGC': '004', 'CGA': '004', 'CGG': '004', 'ATT': '003', 'ATC': '003', 'ATA': '003', 'ATG': '000', 'ACT': '004', 'ACC': '004', 'ACA': '004', 'ACG': '004', 'AAT': '002', 'AAC': '002', 'AAA': '002', 'AAG': '022', 'AGT': '003', 'AGC': '003', 'AGA': '003', 'AGG': '020', 'GTT': '004', 'GTC': '004', 'GTA': '004', 'GTG': '004', 'GCT': '004', 'GCC': '004', 'GCA': '004', 'GCG': '004', 'GAT': '002', 'GAC': '002', 'GAA': '002', 'GAG': '002', 'GGT': '004', 'GGC': '004', 'GGA': '004', 'GGG': '004'}
		CandidateDivisionSR1andGracilibacteria = {

			'TTT': '002', 'TTC': '002', 'TTA': '202', 'TTG': '202', 'TCT': '004', 'TCC': '004', 'TCA': '004', 'TCG': '004', 'TAT': '002', 'TAC': '002', 'TAA': '002', 'TAG': '002', 'TGT': '002', 'TGC': '002', 'TGA': '200', 'TGG': '000', 'CTT': '004', 'CTC': '004', 'CTA': '204', 'CTG': '204', 'CCT': '004', 'CCC': '004', 'CCA': '004', 'CCG': '004', 'CAT': '002', 'CAC': '002', 'CAA': '002', 'CAG': '002', 'CGT': '004', 'CGC': '004', 'CGA': '204', 'CGG': '204', 'ATT': '003', 'ATC': '003', 'ATA': '003', 'ATG': '000', 'ACT': '004', 'ACC': '004', 'ACA': '004', 'ACG': '004', 'AAT': '002', 'AAC': '002', 'AAA': '002', 'AAG': '002', 'AGT': '002', 'AGC': '002', 'AGA': '202', 'AGG': '202', 'GTT': '004', 'GTC': '004', 'GTA': '004', 'GTG': '004', 'GCT': '004', 'GCC': '004', 'GCA': '004', 'GCG': '004', 'GAT': '002', 'GAC': '002', 'GAA': '002', 'GAG': '002', 'GGT': '004', 'GGC': '004', 'GGA': '204', 'GGG': '004'}
		PachysolentannophilusNuclear = {

			'TTT': '002', 'TTC': '002', 'TTA': '202', 'TTG': '002', 'TCT': '004', 'TCC': '004', 'TCA': '004', 'TCG': '004', 'TAT': '002', 'TAC': '002', 'TAA': '022', 'TAG': '002', 'TGT': '002', 'TGC': '002', 'TGA': '020', 'TGG': '000', 'CTT': '003', 'CTC': '003', 'CTA': '203', 'CTG': '000', 'CCT': '004', 'CCC': '004', 'CCA': '004', 'CCG': '004', 'CAT': '002', 'CAC': '002', 'CAA': '002', 'CAG': '002', 'CGT': '004', 'CGC': '004', 'CGA': '204', 'CGG': '204', 'ATT': '003', 'ATC': '003', 'ATA': '003', 'ATG': '000', 'ACT': '004', 'ACC': '004', 'ACA': '004', 'ACG': '004', 'AAT': '002', 'AAC': '002', 'AAA': '002', 'AAG': '002', 'AGT': '002', 'AGC': '002', 'AGA': '202', 'AGG': '202', 'GTT': '004', 'GTC': '004', 'GTA': '004', 'GTG': '004', 'GCT': '004', 'GCC': '004', 'GCA': '004', 'GCG': '004', 'GAT': '002', 'GAC': '002', 'GAA': '002', 'GAG': '002', 'GGT': '004', 'GGC': '004', 'GGA': '004', 'GGG': '004'}
		KaryorelictNuclear = {

			'TTT': '002', 'TTC': '002', 'TTA': '202', 'TTG': '202', 'TCT': '004', 'TCC': '004', 'TCA': '004', 'TCG': '004', 'TAT': '002', 'TAC': '002', 'TAA': '202', 'TAG': '202', 'TGT': '002', 'TGC': '002', 'TGA': '002', 'TGG': '002', 'CTT': '004', 'CTC': '004', 'CTA': '204', 'CTG': '204', 'CCT': '004', 'CCC': '004', 'CCA': '004', 'CCG': '004', 'CAT': '002', 'CAC': '002', 'CAA': '202', 'CAG': '202', 'CGT': '004', 'CGC': '004', 'CGA': '204', 'CGG': '204', 'ATT': '003', 'ATC': '003', 'ATA': '003', 'ATG': '000', 'ACT': '004', 'ACC': '004', 'ACA': '004', 'ACG': '004', 'AAT': '002', 'AAC': '002', 'AAA': '002', 'AAG': '002', 'AGT': '002', 'AGC': '002', 'AGA': '202', 'AGG': '202', 'GTT': '004', 'GTC': '004', 'GTA': '004', 'GTG': '004', 'GCT': '004', 'GCC': '004', 'GCA': '004', 'GCG': '004', 'GAT': '002', 'GAC': '002', 'GAA': '002', 'GAG': '002', 'GGT': '004', 'GGC': '004', 'GGA': '004', 'GGG': '004'}
		CondylostomaNuclear = {

			'TTT': '002', 'TTC': '002', 'TTA': '202', 'TTG': '202', 'TCT': '004', 'TCC': '004', 'TCA': '004', 'TCG': '004', 'TAT': '002', 'TAC': '002', 'TAA': '202', 'TAG': '202', 'TGT': '002', 'TGC': '002', 'TGA': '002', 'TGG': '002', 'CTT': '004', 'CTC': '004', 'CTA': '204', 'CTG': '204', 'CCT': '004', 'CCC': '004', 'CCA': '004', 'CCG': '004', 'CAT': '002', 'CAC': '002', 'CAA': '202', 'CAG': '202', 'CGT': '004', 'CGC': '004', 'CGA': '204', 'CGG': '204', 'ATT': '003', 'ATC': '003', 'ATA': '003', 'ATG': '000', 'ACT': '004', 'ACC': '004', 'ACA': '004', 'ACG': '004', 'AAT': '002', 'AAC': '002', 'AAA': '002', 'AAG': '002', 'AGT': '002', 'AGC': '002', 'AGA': '202', 'AGG': '202', 'GTT': '004', 'GTC': '004', 'GTA': '004', 'GTG': '004', 'GCT': '004', 'GCC': '004', 'GCA': '004', 'GCG': '004', 'GAT': '002', 'GAC': '002', 'GAA': '002', 'GAG': '002', 'GGT': '004', 'GGC': '004', 'GGA': '004', 'GGG': '004'}
		MesodiniumNuclear = {

			'TTT': '002', 'TTC': '002', 'TTA': '202', 'TTG': '202', 'TCT': '004', 'TCC': '004', 'TCA': '004', 'TCG': '004', 'TAT': '004', 'TAC': '004', 'TAA': '004', 'TAG': '004', 'TGT': '002', 'TGC': '002', 'TGA': '000', 'TGG': '000', 'CTT': '004', 'CTC': '004', 'CTA': '204', 'CTG': '204', 'CCT': '004', 'CCC': '004', 'CCA': '004', 'CCG': '004', 'CAT': '002', 'CAC': '002', 'CAA': '002', 'CAG': '002', 'CGT': '004', 'CGC': '004', 'CGA': '204', 'CGG': '204', 'ATT': '003', 'ATC': '003', 'ATA': '003', 'ATG': '000', 'ACT': '004', 'ACC': '004', 'ACA': '004', 'ACG': '004', 'AAT': '002', 'AAC': '002', 'AAA': '002', 'AAG': '002', 'AGT': '002', 'AGC': '002', 'AGA': '202', 'AGG': '202', 'GTT': '004', 'GTC': '004', 'GTA': '004', 'GTG': '004', 'GCT': '004', 'GCC': '004', 'GCA': '004', 'GCG': '004', 'GAT': '002', 'GAC': '002', 'GAA': '002', 'GAG': '002', 'GGT': '004', 'GGC': '004', 'GGA': '004', 'GGG': '004'}
		PeritrichNuclear = {

			'TTT': '002', 'TTC': '002', 'TTA': '202', 'TTG': '202', 'TCT': '004', 'TCC': '004', 'TCA': '004', 'TCG': '004', 'TAT': '002', 'TAC': '002', 'TAA': '202', 'TAG': '202', 'TGT': '002', 'TGC': '002', 'TGA': '000', 'TGG': '000', 'CTT': '004', 'CTC': '004', 'CTA': '204', 'CTG': '204', 'CCT': '004', 'CCC': '004', 'CCA': '004', 'CCG': '004', 'CAT': '002', 'CAC': '002', 'CAA': '002', 'CAG': '002', 'CGT': '004', 'CGC': '004', 'CGA': '204', 'CGG': '204', 'ATT': '003', 'ATC': '003', 'ATA': '003', 'ATG': '000', 'ACT': '004', 'ACC': '004', 'ACA': '004', 'ACG': '004', 'AAT': '002', 'AAC': '002', 'AAA': '002', 'AAG': '002', 'AGT': '002', 'AGC': '002', 'AGA': '202', 'AGG': '202', 'GTT': '004', 'GTC': '004', 'GTA': '004', 'GTG': '004', 'GCT': '004', 'GCC': '004', 'GCA': '004', 'GCG': '004', 'GAT': '002', 'GAC': '002', 'GAA': '202', 'GAG': '202', 'GGT': '004', 'GGC': '004', 'GGA': '004', 'GGG': '004'}
		BlastocrithidiaNuclear = {

			'TTT': '002', 'TTC': '002', 'TTA': '202', 'TTG': '202', 'TCT': '004', 'TCC': '004', 'TCA': '004', 'TCG': '004', 'TAT': '002', 'TAC': '002', 'TAA': '202', 'TAG': '202', 'TGT': '002', 'TGC': '002', 'TGA': '002', 'TGG': '002', 'CTT': '004', 'CTC': '004', 'CTA': '204', 'CTG': '204', 'CCT': '004', 'CCC': '004', 'CCA': '004', 'CCG': '004', 'CAT': '002', 'CAC': '002', 'CAA': '002', 'CAG': '002', 'CGT': '004', 'CGC': '004', 'CGA': '204', 'CGG': '204', 'ATT': '003', 'ATC': '003', 'ATA': '003', 'ATG': '000', 'ACT': '004', 'ACC': '004', 'ACA': '004', 'ACG': '004', 'AAT': '002', 'AAC': '002', 'AAA': '002', 'AAG': '002', 'AGT': '002', 'AGC': '002', 'AGA': '202', 'AGG': '202', 'GTT': '004', 'GTC': '004', 'GTA': '004', 'GTG': '004', 'GCT': '004', 'GCC': '004', 'GCA': '004', 'GCG': '004', 'GAT': '002', 'GAC': '002', 'GAA': '202', 'GAG': '202', 'GGT': '004', 'GGC': '004', 'GGA': '004', 'GGG': '004'}
		BalanophoraceaePlastid = {

			'TTT': '002', 'TTC': '002', 'TTA': '202', 'TTG': '202', 'TCT': '004', 'TCC': '004', 'TCA': '004', 'TCG': '004', 'TAT': '002', 'TAC': '002', 'TAA': '000', 'TAG': '000', 'TGT': '002', 'TGC': '002', 'TGA': '000', 'TGG': '000', 'CTT': '004', 'CTC': '004', 'CTA': '204', 'CTG': '204', 'CCT': '004', 'CCC': '004', 'CCA': '004', 'CCG': '004', 'CAT': '002', 'CAC': '002', 'CAA': '002', 'CAG': '002', 'CGT': '004', 'CGC': '004', 'CGA': '204', 'CGG': '204', 'ATT': '003', 'ATC': '003', 'ATA': '003', 'ATG': '000', 'ACT': '004', 'ACC': '004', 'ACA': '004', 'ACG': '004', 'AAT': '002', 'AAC': '002', 'AAA': '002', 'AAG': '002', 'AGT': '002', 'AGC': '002', 'AGA': '202', 'AGG': '202', 'GTT': '004', 'GTC': '004', 'GTA': '004', 'GTG': '004', 'GCT': '004', 'GCC': '004', 'GCA': '004', 'GCG': '004', 'GAT': '002', 'GAC': '002', 'GAA': '002', 'GAG': '002', 'GGT': '004', 'GGC': '004', 'GGA': '004', 'GGG': '004'}

		if(self.codonTable == 'standard'):
			degenerateCodonTable = standardDict
		elif(self.codonTable == 'VertebrateMitochondrial'):
			degenerateCodonTable = VertebrateMitochondrial
		elif(self.codonTable == 'YeastMitochondrial'):
			degenerateCodonTable = YeastMitochondrial
		elif(self.codonTable == 'MoldMitochondrialProtozoanMitochondrialCoelenterateMitochondrialMycoplasmaSpiroplasma'):
			degenerateCodonTable = MoldMitochondrialProtozoanMitochondrialCoelenterateMitochondrialMycoplasmaSpiroplasma
		elif(self.codonTable == 'InvertebrateMitochondrial'):
			degenerateCodonTable = InvertebrateMitochondrial
		elif(self.codonTable == 'CiliateNuclearDasycladaceanNuclearHexamitaNuclear'):
			degenerateCodonTable = CiliateNuclearDasycladaceanNuclearHexamitaNuclear
		elif(self.codonTable == 'EchinodermMitochondrialFlatwormMitochondrial'):
			degenerateCodonTable = EchinodermMitochondrialFlatwormMitochondrial
		elif(self.codonTable == 'EuplotidNuclear'):
			degenerateCodonTable = EuplotidNuclear
		elif(self.codonTable == 'BacterialArchaealandPlantPlastid'):
			degenerateCodonTable = BacterialArchaealandPlantPlastid
		elif(self.codonTable == 'AlternativeYeastNuclear'):
			degenerateCodonTable = AlternativeYeastNuclear
		elif(self.codonTable == 'AlternativeFlatwormMitochondrial'):
			degenerateCodonTable = AlternativeFlatwormMitochondrial
		elif(self.codonTable == 'BlepharismaMacronuclear'):
			degenerateCodonTable = BlepharismaMacronuclear
		elif(self.codonTable == 'ChlorophyceanMitochondrial'):
			degenerateCodonTable = ChlorophyceanMitochondrial
		elif(self.codonTable == 'TrematodeMitochondrial'):
			degenerateCodonTable = TrematodeMitochondrial
		elif(self.codonTable == 'ScenedesmusobliquusMitochondrial'):
			degenerateCodonTable = ScenedesmusobliquusMitochondrial	
		elif(self.codonTable == 'ThraustochytriumMitochondrial'):
			degenerateCodonTable = ThraustochytriumMitochondrial
		elif(self.codonTable == 'PterobranchiaMitochondrial'):
			degenerateCodonTable = PterobranchiaMitochondrial	
		elif(self.codonTable == 'CandidateDivisionSR1andGracilibacteria'):
			degenerateCodonTable = CandidateDivisionSR1andGracilibacteria
		elif(self.codonTable == 'PachysolentannophilusNuclear'):
			degenerateCodonTable = PachysolentannophilusNuclear
		elif(self.codonTable == 'KaryorelictNuclear'):
			degenerateCodonTable = KaryorelictNuclear
		elif(self.codonTable == 'CondylostomaNuclear'):
			degenerateCodonTable = CondylostomaNuclear
		elif(self.codonTable == 'MesodiniumNuclear'):
			degenerateCodonTable = MesodiniumNuclear
		elif(self.codonTable == 'PeritrichNuclear'):
			degenerateCodonTable = PeritrichNuclear
		elif(self.codonTable == 'BlastocrithidiaNuclear'):
			degenerateCodonTable = BlastocrithidiaNuclear
		elif(self.codonTable == 'BalanophoraceaePlastid'):
			degenerateCodonTable = BalanophoraceaePlastid



		degenerancy = ''
		for i in range(0, len(data),3):
			codon = data[i:i+3]
			if('N' in codon or '-' in codon):
				degenerancy += codon
			else:
				degenerancy += degenerateCodonTable[codon]

		return(degenerancy)
	
	def sequencesToMatrix(self):

		fasta = px.Fasta(self.file,duplicate_action='first',sequence_always_upper=True,read_long_names=True)

		# Extract samples from fastas
		samples = list(fasta.keys())
		
		seqLen = len(fasta[samples[0]][:].seq)
		
		if((seqLen % 3) != 0):
			print('cdsLength')
			sys.exit('cdsLength')

		# Create empty array with ndimesions equal to multi-Fasta lines and length
		matrix = np.empty([len(samples),len(fasta[samples[0]][:].seq)],dtype='str')
		
		# List to append indexes if whole sequence at any population is len(seq) * 'N'
		deleteIndex = list()
		
		# Iter fasta to add sequence to matrix
		for i in range(1,len(samples),1):
			# Extract each sample sequence
			tmp = fasta[samples[i]][:].seq
			if(len(tmp) != seqLen):
				print('errorAlign')
				sys.exit('errorAlign')

			# if('N' in tmp):
			if('N'*len(tmp) == tmp):
				deleteIndex.append(i)
			else:
				matrix[i] = list(tmp)
		
		# Delete lines
		matrix = np.delete(matrix,deleteIndex,0)
		
		degenCode = self.degenerancy(fasta[samples[0]][:].seq)
		# Put degenerancy in first ndarray element
		matrix[0] = list(degenCode)
		# NEED TO SOLVE THIS. C-contigous change web subset, need true in order to inter properly witih nditer

		matrix = np.asarray(matrix[:,(matrix[0]=='0') | (matrix[0]=='4')],order='C')

		return(matrix)

	def uSfsFromFasta(self, matrix):

		output = list()
		for x in np.nditer(matrix, order='F',flags=['external_loop']): 
			degen = x[0]
			AA = x[-1]

			# Undefined Ancestra Allele. Try to clean out of the loop
			if(AA == 'N' or AA == '-'):
				next
			elif('N' in x[1:-1] or '-' in x[1:-1]):
				next
			# Monomorphic sites. Try to clean out of the loop
			elif(np.unique(x[1:][np.where(x[1:]!='N')]).shape[0] == 1):
				next
			else:
				pol = x[1:-1]
				if(degen == '4'):
					functionalClass = '4fold'
				else:
					functionalClass = '0fold'

				# Check if pol != AA and monomorphic
				if((np.unique(pol).shape[0] == 1) and (np.unique(pol)[0] != AA)):
					div = 1; AF = 0
					tmp = [AF,div,functionalClass]
					output.append(tmp)
				else:
					AN = x[1:-1].shape[0]
					AC = pd.DataFrame(data=np.unique(x[1:-1], return_counts=True)[1],index=np.unique(x[1:-1], return_counts=True)[0])
					div = 0
					if(AA not in AC.index):
						next
					else:
						AC = AC[AC.index!=AA]
						if(len(AC) == 0):
							next
						else:
							AF = AC.iloc[0]/AN
							AF = AF.iloc[0]
					tmp = [AF,div,functionalClass]
					output.append(tmp)
		return(output)

	def formatSfs(self, sfs, matrix):

		df         = pd.DataFrame(sfs)
		df['id']   = 'uploaded'
		df.columns = ['derivedAlleleFrequency','d','functionalClass','id']

		# Extract divergence data
		div = df[['id','functionalClass','d']]
		div = div[div['d']!=0]
		div = div.groupby(['id','functionalClass'])['d'].count().reset_index()
		div = div.pivot_table(index=['id'],columns=['functionalClass'],values='d').reset_index()
		try:
			div = div[['0fold','4fold']]
		except:
			if('4fold' in div.columns):
				div = div[['4fold']]
				div['0fold'] = 0
			elif('0fold' in div.columns):
				div = div[['0fold']]
				div['4fold'] = 0

		div['mi'] =  matrix[0][np.where(matrix[0]=='0')].shape[0]
		div['m0'] =  matrix[0][np.where(matrix[0]=='4')].shape[0]
		div = div.rename(columns={'0fold':'Di','4fold':'D0','mi':'mi','m0':'m0'})
		# div = div.pivot_table(index=['functionalClass'],columns=['functionalClass'],values='div').reset_index()

		# Create SFS pd.DataFrame by functionClass and 20 frequency bin
		daf = df[df['d']!=1][['derivedAlleleFrequency','functionalClass','id']]
		b = np.linspace(0,1,self.bins)
		labels = b[1:].tolist()
		daf['categories'] = pd.cut(daf['derivedAlleleFrequency'],bins=b,labels=labels)
		daf = daf.groupby(['functionalClass','id','categories']).count().reset_index()
		sfs = pd.DataFrame({'daf':daf['categories'].unique(),'P0':daf[daf['functionalClass']=='4fold']['derivedAlleleFrequency'].reset_index(drop=True),'Pi':daf[daf['functionalClass']=='0fold']['derivedAlleleFrequency'].reset_index(drop=True)})

		sfs = sfs[['daf','P0','Pi']]
		sfs['P0'] = sfs['P0'].fillna(0)
		sfs['Pi'] = sfs['Pi'].fillna(0)
		sfs['daf'] = sfs['daf'].apply(lambda x: round(x,3))

		return(np.array(sfs),np.array(div))