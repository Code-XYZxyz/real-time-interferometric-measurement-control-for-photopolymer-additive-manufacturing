---
title: 'A MATLAB-based software for realizing an affordable real-time measured and controlled ultraviolet curing system in photopolymer additive manufacturing'
tags:
  - Additive manufacturing
  - Photopolymerization
  - Real-time feedback control
  - Process measurement
authors:
 - name: Xiayun Zhao
   orcid: 0000-0002-6092-1130
   affiliation: 1
 - name: David W. Rosen
   orcid: N/A
   affiliation: 1
affiliations:
 - name: Georgia Institute of Technology
   index: 1
date: 12 September 2017
bibliography: paper.bib
---

# Summary

- This is a comprehensive MATLAB-based software platform developed for real-time measurement and feedback control of a custom mask-projection photopolymerization based additive manufacturing system (referred as "ECPL", i.e., Exposure Controlled Projection Lithography) using a lab-built interferometry (referred as "ICM&M", i.e., Interferometric Curing Monitoring and Measurement). A graphical user interface using the graphical user interface development environment (GUIDE) of MATLAB was created to implement the ICM&M method for the ECPL process. The software interfaces with the hardware of the ECPL system’s ultraviolet lamp and DMD, and the ICM&M system’s camera. It was designed to streamline the operation of the ECPL process with the aid of parallel computing that implements online both the ICM&M acquisition and measurement analysis as well as the feedback control method. The application logs the acquired interferogram video data, performs numerical computations for the ICM&M measurement algorithms and control law, saves the real-time data and measurement results for all voxels in the region of interest. Meanwhile, it displays interferogram frames and visualize the photocuring process without a substantial sacrifice in temporal performance of other key functions such as data acquisition and measurement & control analysis. 

- The software could be extended to real-time process measurement and control for other additive manufacturing systems, for example, stereo-lithography (SLA) and metal based additive manufacturing aided by in-situ thermal images analysis.

- Related reference publications about the ECPL additive manufacturing system design, the ICM&M principle and validation, the real-time experiment results are listed below (in the references section). The software presented here is the backbone of the physical implementation of the ECPL process measurement and feedback control.

- Demo videos of how to use the software, examples and corresponding metadata are provided in this repository’s folder of “Examples-Metadata”.

- The document “paper.pdf” contains: (1) the research application that is associated with the software; (2) details about the software design, functions, and flowchart; (3) implementation examples.


# References

1.	X. Zhao and D. Rosen, “Real-time Interferometric Monitoring and Measuring of photopolymerization based stereolithographic additive manufacturing process: sensor model and algorithm”, Measurement Science and Technology, Vol 28, Issue 1, 2016. doi:10.1088/0957-0233/28/1/015001. (http://dx.doi.org/10.1088/0957-0233/28/1/015001).
2.	X. Zhao and D. Rosen, “A data mining approach of real-time process measurement for polymer additive manufacturing with the exposure controlled projection lithography”, Journal of Manufacturing Systems, Special Issue: Cybermanufacturing. doi:10.1016/j.jmsy.2017.01.005. (http://dx.doi.org/10.1016/j.jmsy.2017.01.005).
3.	X. Zhao and D. Rosen, “Experimental validation and characterization of a real-time metrology system for photopolymerization based stereolithographic additive manufacturing process”, International Journal of Advanced Manufacturing Technology, Dec 2016. doi:10.1007/s00170-016-9844-1. (http://dx.doi.org/10.1007/s00170-016-9844-1)
4.	X. Zhao and D. Rosen, “Parallel computing enabled real-time interferometric measurement and feedback control for a photopolymer based lithographic additive manufacturing process”, Mechatronics, Special Issue: Mechatronics and Additive Manufacturing. (July 2017: in review)
5.	X. Zhao, “Process Measurement and Control for Exposure Controlled Projection Lithography”. Ph.D. Dissertation, Mechanical Engineering, Georgia Institute of Technology, Atlanta, USA, 2017. Available on https://smartech.gatech.edu/handle/1853/58294
6.	X. Zhao and D. Rosen, “Simulation Study on Evolutionary Cycle to Cycle Time Control of Exposure Controlled Projection Lithography”, Rapid Prototyping Journal, Vol 22, Issue 3, 2016, pp 456-464. doi: 10.1108/RPJ-01-2015-0008. doi:10.1108/RPJ-01-2015-0008. (http://dx.doi.org/10.1108/RPJ-01-2015-0008)
