# EigenComputing
SVD problem solving and eigenvalue calculation using the Jacobi method.

Two Jacobi methods were implemented in the Julia programming language to obtain eigenvalues and eigenvectors of symmetric matrices. One is called -Jacobi Classic-, and the other is called -Cyclic by row-. One is called  and the other is called. An examples are shown below,

For Jacobi Classic
![Picture1](https://github.com/YairCCastillo/EigenComputing/assets/49602985/84bf7c60-9a49-4590-875d-a863796d09a7)

For Cyclic by row
![Picture2](https://github.com/YairCCastillo/EigenComputing/assets/49602985/c3376b47-f1f8-42ee-8e1e-c0fe6a02823d)

where D is the matrix of eigenvalues on the diagonal and B is the matrix of eigenvectors.

Also we compare the average number of sweeps to reach the eigenvalus and eigenvectors:
![image](https://github.com/YairCCastillo/EigenComputing/assets/49602985/d1f3e22e-be4e-49a2-b340-2d13f763ff8f)

## Solving SVD problemas

Something astonishing about the Jacobi method is that we can also extend it to solve SVD problems. 

First example:
![image](https://github.com/YairCCastillo/EigenComputing/assets/49602985/992bda1b-dd09-45fd-85f7-c30c6fb8def7)

Second example:
![image](https://github.com/YairCCastillo/EigenComputing/assets/49602985/d012ce8b-315a-4c05-9836-2ad80ccb5426)
