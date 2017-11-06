 <h1>Implementation of block matrix multiplication</h1>
 

<p> This is impementation of block matrix multiplication consists from 3 types</p>

<ul>
 <li>Sequential block multiplication</li>
 <li>Internal parallelization of two blocks</li>
 <li>External parallelization of sycle of multiplication(When two different blocks multiply on two different cores)</li>
</ul>


<p style="color: blue; font-size: 16pt;"><b>Requirements:</b></p>
<ol>
 <li>Matrix A must be a symmetrical and stored as a low-triangle matrix by rows</li>
 <li>Matrix B must be a top-triangular and stored without zeros by columns</li>
</ol>

<p>Internal and external parallelization based on OpenMP technology</p> 
<br>

<p> This repository contains the results of the tests for GCC and Clang compilers for two different types: double and float, respectively</p>


----------------------------------------------------------------------------------------


<h3> Results of Clang </h3>

![Float clang](https://rawgit.com/rekrut1993/Block-Matrix-Multiplication/master/results/clang/float/clang_float.svg)

![Double clang](https://rawgit.com/rekrut1993/Block-Matrix-Multiplication/master/results/clang/double/clang_double.svg) 


<h3> Results of GCC </h3>

![Float gcc](https://rawgit.com/rekrut1993/Block-Matrix-Multiplication/master/results/gcc/float/gcc_float.svg)

![Double gcc](https://rawgit.com/rekrut1993/Block-Matrix-Multiplication/master/results/gcc/double/gcc_double.svg) 


----------------------------------------------------------------------------------------

<h3> How build </h3>
<ol>
 <li>Clone this repo</li>
 <li>make <i>cd</i> in this repo</li>
 <li><i>mkdir</i> build && <i>cd</i> build/</li>
 <li><i>cmake</i> ../</li>
 <li><i>make</i> && <i>cd</i> ..</li>
 <li>.build/MultiplicationBlockMatrix</li>
</ol>
