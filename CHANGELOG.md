# CHANGELOG
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

### Somatic_whole_genome_sequencing_pipeline-v3.0.0
	- Features
		- Added a new dsl2 structure similar to Somatic Panel Pipeline
		- Added the version information for all the process and softwares 
		- Added the common workflow entry point "SWGP" in order to aactivate the somtic pipeline workflow
		- Seperated the "Hema" and "solid" profile based on the diagnosis question
	- Fixes
		- Added the Structural variants and tandem SV analysis in gene fusion identification
		- Fixed the new gens link for coyote


### Release date 
	* version 2.0.0 on September 30 2022

#### Features
	*	DSL2 based implemetation of tumor WGS
	*	Able to adabpt to both Hematology and Solid WGS questions
	* 	GENS for both normal and tumor
	* 	Able to handle both tumor/normal and tumor only analysis
	* 	COYOTE can handle unfiltered SNV list but tagged SNV from gene panel can be selected
	
	
#### Fixes
	* 	Sharding from 32 to 8 for better I/O specially with bqsr steps in senteion
	* 	Changes with the profile configs
	* 	Changes with the coyote to include sample-wgs
