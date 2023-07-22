COMPUTABLE PERFORMANCE GUARANTEES FOR COMPRESSED SENSING MATRICES
(Subtitle: How to check the recovery performance of a given sensing matrix in compressed sensing?)
------------------------
A. INTRODUCE
The null space condition for L1 minimization in compressed sensing is a necessary and sufficient condition on the sensing matrices under which a sparse signal can be uniquely recovered from the observation data via L1 minimization. However, verifying the null space condition is known to be computationally challenging. Most of the existing methods can only provide upper and lower bounds on the proportion parameter that characterizes the null space condition. In this paper, we propose new polynomial-time algorithms to establish the upper bounds of the proportion parameter. Based on these polynomial-time algorithms, we have designed new algorithms - Sandwiching Algorithm (SWA) and Tree Search Algorithm (TSA) - to precisely verify the null space condition. Extensive simulations show that Tree Search Algorithm and Sandwiching Algorithm significantly reduce the computational complexity when compared with the exhaustive search method.

------------------------
B. GOALS:
Providing the maximum recoverable sparsity k via L1 minimization given a sensing matrix A(size: m x n)

------------------------
C. CONTENTS:
matlab code (MOSEK solver should be installed. (available at https://www.mosek.com/))
!!!! Prerequisite software (or solver) to run the MATLAB code: Mosek, CVX !!!!
https://www.mosek.com
http://cvxr.com


C.1. FOLDERS
00_Rawdata: raw data file for sensing matrices used to use for Table 1 and 2 of [1].
01_JNcode: LP relaxation code obtained from http://www2.isye.gatech.edu/~nemirovs/
02_Pick-l: pick-l algorithm code
03_SWA_StandAlone: SandWiching Algorithm(SWA) (introduced in sconference processing [2])
04_TSA_StandAlone: Tree Search Algorithm(TSA)

C.2. MAIN FILES
T0_NSP_SWA.m" Testing Sandwiching Algorithm
T2_NSP_HD_comp.m: lower bound on recoverable sparsity k on high-dimensional sensing matrix (n=1024,2048,etc.)(Table 2 of [1])
T4_NSP_LD_comp_Gauss_optPick: Regenerating Table 4 of [1]
T5_NSP_LD_comp_Gauss.m: upper bound on alpha_k comparision on low-dimensional Gaussian sensing matrix (n=40)(Table 5 of [1])
T6_NSP_LD_comp_Four.m: upper bound on alpha_k comparision on low-dimensional Fourier sensing matrix (n=40)(Table 5 of [1])
T7_NSP_LD_comp_Gauss.m: Regenerating Table 7 of [1]
T8_NSP_LD_comp_Four.m: Regenerating Table 8 of [1]
F5_NSP_TSA_Complexity_Gauss.m: Regenerating Figure 5 of [1] for the complexity of TSA on Gaussian matrix


------------------------
D. ADDITION
(a) To be consistent with the previous research, we used the Matlab routine [2] provided at http://www2.isye.gatech.edu/nemirovs/ for generating sensing matrices. (AAGenerateSensingMatrix.m)

(b) Contact info.: Myung (Michael) Cho 
- Email: michael.myung.cho@gmail.com
- Homepage: https://sites.google.com/view/myungcho/home

------------------------
E. REFERENCE PAPER: 
[1] M. Cho, K. V. Mishra, and W. Xu, "Computable performance guarantees for compressed sensing matrices," EURASIP Journal on Advances in Signal Processing, no. 16, 2018.
https://asp-eurasipjournals.springeropen.com/articles/10.1186/s13634-018-0535-y

[2] M. Cho, W. Xu, "New algorithms for verifying the null space conditions in compressed sensing," In Proceedings of
Asilomar Conference on Signals, Systems and Computers (2013), pp. 1038-1042.
https://ieeexplore.ieee.org/document/6810449

[3] A. Juditsky and A. Nemirovski, “On verifiable sufficient conditions for sparse signal recovery via L1 minimization,” Mathematical programming, vol. 127, no. 1, pp. 57–88, 2011.
[4] A. d’Aspremont and L. El Ghaoui, “Testing the nullspace property using semidefinite programming,” Mathematical programming, vol. 127, no. 1, pp. 123–144, 2011.


------------------------
F. FAQ:
1. If "Undefined function or variable 'find_portion'" error occured, then run the following commend to build mex file.
(matlab prompt)>>./04_TSA_StandAlone/mex find_portion.c
