import { sum } from 'mathjs';

// Declaring global variables
let Nx = 50;
let Ny = 50;
let Re = 100;
let Mew = 0.1;
let tStop = 10;
    
let u = Array.from({length: Nx}, () => 
    Array(Ny).fill(1)
);
let v = Array.from({length: Nx}, () => 
    Array(Ny).fill(1)
);
let rho = Array.from({length: Nx}, () => 
    Array(Ny).fill(1)
);

const Ksi = [
    [0, 1, 0, -1, 0, 1, -1, -1, 1],
    [0, 0, 1, 0, -1, 1, 1, -1, -1]
];

const w = [4/9, 1/9, 1/9, 1/9, 1/9, 1/36, 1/36, 1/36/ 1/36];

const cs = 1/(Math.sqrt(3));

const dx = 1;
const dy = 1;
const L = (Nx - 1)*dx;

const uLid = Re*Mew/L;

const Tau = 3*Mew+0.5;

let fNew = Array.from({length: Ny}, () => 
    Array.from({ length: Nx }, () => 
        Array(9).fill(1)
    )
);
let fOld = Array.from({length: Ny}, () => 
    Array.from({ length: Nx }, () => 
        Array(9).fill(1)
    )
);

let fEq = Array.from({length: Ny}, () => 
    Array.from({ length: Nx }, () => 
        Array(9).fill(1)
    )
);

let rhoB;

function calculate(){
    // Calculation
    for (let t = 0; t < tStop; t++){
        // Streaming for all interior nodes and for boundary nodes
        for (let m = 0; m < Nx; m++){ // x coordinates
            for (let n = 0; n < Ny; n++){ // y coordinates
                if (n == 0){ // top boundary
                    if (m == 0){
                        topLeftCornerNode(n,m);
                    }else if (m == (Nx - 1)){
                    topRightCornerNode(n,m);
                }else{
                    otherTopNodes(n,m);
                }
            }else if(n == (Ny - 1)){ // bottom boundary
                if (m == 0){
                    bottomLeftCornerNode(n,m);
                }else if (m == (Nx - 1)){
                    bottomRightCornerNode(n,m);
                }else{
                    otherBottomNodes(n,m);
                }
            }else if (m == 0){ // left boundary
                leftBoundary(n,m);
            }else if (m == (Nx - 1)){ // right boundary
                rightBoundary(n,m);
            }else{
                interiorNodes(n,m);
            }
            
        }
        for (let m = 0; m < Nx; m++){
            for (let n = 0; n < Ny; n++){
                momentCalculation(n,m);
            }
        }
        for (let m = 0; m < Nx; m++){
            for (let n = 0; n < Ny; n++){
                fEqCalculation(n,m);
            }
        }
    }
}
}

function topLeftCornerNode(n,m){
    //console.log("left");
    fNew[n,m,0] = fOld[n,m,0];
    fNew[n,m,2] = fOld[n+1,m,2];
    fNew[n,m,3] = fOld[n,m+1,3];
    fNew[n,m,6] = fOld[n+1,m+1,6];

    fNew[n,m,1] = fNew[n,m,3];
    fNew[n,m,4] = fNew[n,m,2];
    fNew[n,m,8] = fNew[n,m,6];

    rhoB = (rho[n+1,m] + rho[n,m+1])/2;
    fNew[n,m,5] = (rhoB-fNew[n,m,0] - fNew[n,m,1]-fNew[n,m,2]-fNew[n,m,3]-fNew[n,m,4]-fNew[n,m,6]-fNew[n,m,8])/2;
    fNew[n,m,7] = fNew[n,m,5];
}

function topRightCornerNode(n,m){
    //console.log("right")
    fNew[n,m,0] = fOld[n,m,0];
    fNew[n,m,1] = fOld[n,m-1,1];
    fNew[n,m,2] = fOld[n+1,m,2];
    fNew[n,m,5] = fOld[n+1,m-1,5];

    fNew[n,m,3] = fNew[n,m,1];
    fNew[n,m,4] = fNew[n,m,2];
    fNew[n,m,7] = fNew[n,m,5];

    rhoB = (rho[n+1,m]+rho[n,m-1])/2;
    fNew[n,m,6] = (rhoB-fNew[n,m,0] - fNew[n,m,1] - fNew[n,m,2] - fNew[n,m,3] - fNew[n,m,4] - fNew[n,m,5] - fNew[n,m,7])/2;
    fNew[n,m,8] = fNew[n,m,6];
}

function otherTopNodes(n,m){
    //console.log("other");
    fNew[n,m,0] = fOld[n,m,0];
    fNew[n,m,1] = fOld[n,m-1,1];
    fNew[n,m,2] = fOld[n+1,m,2];
    fNew[n,m,3] = fOld[n,m+1,3];
    fNew[n,m,5] = fOld[n+1,m-1,5];
    fNew[n,m,6] = fOld[n+1,m+1,6];

    fNew[n,m,4] = fNew[n,m,2];
    rhoB = fNew[n,m,0] + fNew[n,m,1] + fNew[n,m,3] + 2*(fNew[n,m,2] + fNew[n,m,5] + fNew[n,m,6]);
    fNew[n,m,7] = fNew[n,m,5] + (fNew[n,m,1] - fNew[n,m,3])/2 - rhoB*uLid/2;
    fNew[n,m,8] = fNew[n,m,6] + (fNew[n,m,3] - fNew[n,m,1])/2 + rhoB*uLid/2;
}

function bottomLeftCornerNode(n,m){
    fNew[n,m,0] = fOld[n,m,0];
    fNew[n,m,3] = fOld[n,m+1,3];
    fNew[n,m,4] = fOld[n-1,m,4];
    fNew[n,m,7] = fOld[n-1,m+1,7];

    fNew[n,m,1] = fNew[n,m,3];
    fNew[n,m,2] = fNew[n,m,4];
    fNew[n,m,5] = fNew[n,m,7];
    rhoB = (rho[n-1,m] + rho[n,m+1])/2;
    fNew[n,m,6] = (rhoB - fNew[n,m,0] - fNew[n,m,1] - fNew[n,m,2] - fNew[n,m,3] - fNew[n,m,4] - fNew[n,m,5] - fNew[n,m,7])/2;
    fNew[n,m,8] = fNew[n,m,6];
}

function bottomRightCornerNode(n,m){
    fNew[n,m,0] = fOld[n,m,0];
    fNew[n,m,1] = fOld[n,m-1,1];
    fNew[n,m,4] = fOld[n-1,m,4];
    fNew[n,m,8] = fOld[n-1,m-1,8];

    fNew[n,m,3] = fNew[n,m,1];
    fNew[n,m,2] = fNew[n,m,4];
    fNew[n,m,6] = fNew[n,m,8];
    rhoB = (rho[n-1,m] + rho[n,m-1])/2;
    fNew[n,m,5] = (rhoB - fNew[n,m,0] - fNew[n,m,1] - fNew[n,m,2] - fNew[n,m,3] - fNew[n,m,4] - fNew[n,m,6] - fNew[n,m,8])/2;
    fNew[n,m,7] =  fNew[n,m,5];
} 

function otherBottomNodes(n,m){
    fNew[n,m,0] = fOld[n,m,0];
    fNew[n,m,1] = fOld[n,m-1,1];
    fNew[n,m,3] = fOld[n,m+1,3];
    fNew[n,m,4] = fOld[n-1,m,4];
    fNew[n,m,7] = fOld[n-1,m+1,7];
    fNew[n,m,8] = fOld[n-1,m-1,8];

    fNew[n,m,2] = fNew[n,m,4];
    fNew[n,m,5] = fNew[n,m,7] + (fNew[n,m,3] - fNew[n,m,1])/2;
    fNew[n,m,6] = fNew[n,m,8] + (fNew[n,m,1] - fNew[n,m,3])/2;
}

function leftBoundary(n,m){
    fNew[n,m,0] = fOld[n,m,0];
    fNew[n,m,2] = fOld[n+1,m,2];
    fNew[n,m,3] = fOld[n,m+1,3];
    fNew[n,m,4] = fOld[n-1,m,4];
    fNew[n,m,6] = fOld[n+1,m+1,6];
    fNew[n,m,7] = fOld[n-1,m+1,7];

    fNew[n,m,1] = fNew[n,m,3];
    fNew[n,m,5] = fNew[n,m,7] + (fNew[n,m,4] - fNew[n,m,2])/2;
    fNew[n,m,8] = fNew[n,m,6] + (fNew[n,m,2] - fNew[n,m,4])/2;
}

function rightBoundary(n,m){
    fNew[n,m,0] = fOld[n,m,0]; 
    fNew[n,m,1] = fOld[n,m-1,1];
    fNew[n,m,2] = fOld[n+1,m,2];
    fNew[n,m,4] = fOld[n-1,m,4];
    fNew[n,m,5] = fOld[n+1,m-1,5];
    fNew[n,m,8] = fOld[n-1,m-1,8];

    fNew[n,m,3] = fNew[n,m,1];
    fNew[n,m,6] = fNew[n,m,8] + (fNew[n,m,4] - fNew[n,m,2])/2;
    fNew[n,m,7] = fNew[n,m,5] + (fNew[n,m,2] - fNew[n,m,4])/2;
}

function interiorNodes(n,m){
    fNew[n,m,0] = fOld[n,m,0];
    fNew[n,m,1] = fOld[n,m-1,1];
    fNew[n,m,2] = fOld[n+1,m,2];
    fNew[n,m,3] = fOld[n,m+1,3];
    fNew[n,m,4] = fOld[n-1,m,4];
    fNew[n,m,5] = fOld[n+1,m-1,5];
    fNew[n,m,6] = fOld[n+1,m+1,6];
    fNew[n,m,7] = fOld[n-1,m+1,7];
    fNew[n,m,8] = fOld[n-1,m-1,8];
}

function momentCalculation(n,m){
    //console.log("m: " + m + ", n: " + n);
    //console.log(fNew[n][m]);
    rho[n][m] = sum(fNew[n][m].slice(0, 9));
    u[n,m] = (fNew[n,m,1] + fNew[n,m,5] + fNew[n,m,8] - fNew[n,m,3] - fNew[n,m,6] - fNew[n,m,7])/rho[n,m];
    v[n,m] = (fNew[n,m,2] + fNew[n,m,5] + fNew[n,m,6] - fNew[n,m,4] - fNew[n,m,7] - fNew[n,m,8])/rho[n,m];
}

function fEqCalculation(n,m){

}

//calculate();

for (let m = 0; m < Nx; m++){
    for (let n = 0; n < Ny; n++){
        console.log(fNew[n,m]);
    }
}