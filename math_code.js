//import { sum } from 'mathjs';

// Declaring global variables
let scaleCell = 10; // Size of each cell on screen
let scaleColor = 1000 // Controls saturation of color mapping
let drawStep = 1; // How often to draw frames
let stepCount = 0;

let Nx = 35;
let Ny = 35;
let Re = 100;
let Mew = 0.1;
let tStop = 20;
    
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
  [0, 0],
  [1, 0],
  [0, 1],
  [-1, 0],
  [0, -1],
  [1, 1],
  [-1, 1],
  [-1, -1],
  [1, -1]
];

const w = [4/9, 1/9, 1/9, 1/9, 1/9, 1/36, 1/36, 1/36, 1/36];

const cs = 1/(Math.sqrt(3));

const dx = 1;
const dy = 1;
const L = (Nx - 1)*dx;

const uLid = Re*Mew/L;

const Tau = 3*Mew+0.5;

let fNew = Array.from({length: Nx}, () => 
    Array.from({ length: Ny }, () => 
        Array(9).fill(1)
    )
);
let fOld = Array.from({length: Nx}, () => 
    Array.from({ length: Ny }, () => 
        Array(9).fill(1)
    )
);

let fEq = Array.from({length: Nx}, () => 
    Array.from({ length: Ny }, () => 
        Array(9).fill(1)
    )
);

let rhoB;

// Calculating Node Data
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
        //console.log("t: " + t);
        for (let m = 0; m < Nx; m++){
            for (let n = 0; n < Ny; n++){
                for(let k = 0; k < 9; k++){
                    fEqCalculation(n,m,k);
                }
            }
        }

        collision();
    }
}
}

function topLeftCornerNode(n,m){
    fNew[n][m][0] = fOld[n][m][0];
    fNew[n][m][2] = fOld[n+1][m][2];
    fNew[n][m][3] = fOld[n][m+1][3];
    fNew[n][m][6] = fOld[n+1][m+1][6];

    fNew[n][m][1] = fNew[n][m][3];
    fNew[n][m][4] = fNew[n][m][2];
    fNew[n][m][8] = fNew[n][m][6];

    rhoB = (rho[n+1][m] + rho[n][m+1])/2;
    fNew[n][m][5] = (rhoB-fNew[n][m][0] - fNew[n][m][1]-fNew[n][m][2]-fNew[n][m][3]-fNew[n][m][4]-fNew[n][m][6]-fNew[n][m][8])/2;
    fNew[n][m][7] = fNew[n][m][5];
}

function topRightCornerNode(n,m){
    fNew[n][m][0] = fOld[n][m][0];
    fNew[n][m][1] = fOld[n][m-1][1];
    fNew[n][m][2] = fOld[n+1][m][2];
    fNew[n][m][5] = fOld[n+1][m-1][5];

    fNew[n][m][3] = fNew[n][m][1];
    fNew[n][m][4] = fNew[n][m][2];
    fNew[n][m][7] = fNew[n][m][5];

    rhoB = (rho[n+1][m]+rho[n][m-1])/2;
    fNew[n][m][6] = (rhoB-fNew[n][m][0] - fNew[n][m][1] - fNew[n][m][2] - fNew[n][m][3] - fNew[n][m][4] - fNew[n][m][5] - fNew[n][m][7])/2;
    fNew[n][m][8] = fNew[n][m][6];
}

function otherTopNodes(n,m){
    fNew[n][m][0] = fOld[n][m][0];
    fNew[n][m][1] = fOld[n][m-1][1];
    fNew[n][m][2] = fOld[n+1][m][2];
    fNew[n][m][3] = fOld[n][m+1][3];
    fNew[n][m][5] = fOld[n+1][m-1][5];
    fNew[n][m][6] = fOld[n+1][m+1][6];

    fNew[n][m][4] = fNew[n][m][2];
    rhoB = fNew[n][m][0] + fNew[n][m][1] + fNew[n][m][3] + 2*(fNew[n][m][2] + fNew[n][m][5] + fNew[n][m][6]);
    fNew[n][m][7] = fNew[n][m][5] + (fNew[n][m][1] - fNew[n][m][3])/2 - rhoB*uLid/2;
    fNew[n][m][8] = fNew[n][m][6] + (fNew[n][m][3] - fNew[n][m][1])/2 + rhoB*uLid/2;
}

function bottomLeftCornerNode(n,m){
    fNew[n][m][0] = fOld[n][m][0];
    fNew[n][m][3] = fOld[n][m+1][3];
    fNew[n][m][4] = fOld[n-1][m][4];
    fNew[n][m][7] = fOld[n-1][m+1][7];

    fNew[n][m][1] = fNew[n][m][3];
    fNew[n][m][2] = fNew[n][m][4];
    fNew[n][m][5] = fNew[n][m][7];
    rhoB = (rho[n-1][m] + rho[n][m+1])/2;
    fNew[n][m][6] = (rhoB - fNew[n][m][0] - fNew[n][m][1] - fNew[n][m][2] - fNew[n][m][3] - fNew[n][m][4] - fNew[n][m][5] - fNew[n][m][7])/2;
    fNew[n][m][8] = fNew[n][m][6];
}

function bottomRightCornerNode(n,m){
    fNew[n][m][0] = fOld[n][m][0];
    fNew[n][m][1] = fOld[n][m-1][1];
    fNew[n][m][4] = fOld[n-1][m][4];
    fNew[n][m][8] = fOld[n-1][m-1][8];

    fNew[n][m][3] = fNew[n][m][1];
    fNew[n][m][2] = fNew[n][m][4];
    fNew[n][m][6] = fNew[n][m][8];
    rhoB = (rho[n-1][m] + rho[n][m-1])/2;
    fNew[n][m][5] = (rhoB - fNew[n][m][0] - fNew[n][m][1] - fNew[n][m][2] - fNew[n][m][3] - fNew[n][m][4] - fNew[n][m][6] - fNew[n][m][8])/2;
    fNew[n][m][7] =  fNew[n][m][5];
} 

function otherBottomNodes(n,m){
    fNew[n][m][0] = fOld[n][m][0];
    fNew[n][m][1] = fOld[n][m-1][1];
    fNew[n][m][3] = fOld[n][m+1][3];
    fNew[n][m][4] = fOld[n-1][m][4];
    fNew[n][m][7] = fOld[n-1][m+1][7];
    fNew[n][m][8] = fOld[n-1][m-1][8];

    fNew[n][m][2] = fNew[n][m][4];
    fNew[n][m][5] = fNew[n][m][7] + (fNew[n][m][3] - fNew[n][m][1])/2;
    fNew[n][m][6] = fNew[n][m][8] + (fNew[n][m][1] - fNew[n][m][3])/2;
}

function leftBoundary(n,m){
    fNew[n][m][0] = fOld[n][m][0];
    fNew[n][m][2] = fOld[n+1][m][2];
    fNew[n][m][3] = fOld[n][m+1][3];
    fNew[n][m][4] = fOld[n-1][m][4];
    fNew[n][m][6] = fOld[n+1][m+1][6];
    fNew[n][m][7] = fOld[n-1][m+1][7];

    fNew[n][m][1] = fNew[n][m][3];
    fNew[n][m][5] = fNew[n][m][7] + (fNew[n][m][4] - fNew[n][m][2])/2;
    fNew[n][m][8] = fNew[n][m][6] + (fNew[n][m][2] - fNew[n][m][4])/2;
}

function rightBoundary(n,m){
    fNew[n][m][0] = fOld[n][m][0]; 
    fNew[n][m][1] = fOld[n][m-1][1];
    fNew[n][m][2] = fOld[n+1][m][2];
    fNew[n][m][4] = fOld[n-1][m][4];
    fNew[n][m][5] = fOld[n+1][m-1][5];
    fNew[n][m][8] = fOld[n-1][m-1][8];

    fNew[n][m][3] = fNew[n][m][1];
    fNew[n][m][6] = fNew[n][m][8] + (fNew[n][m][4] - fNew[n][m][2])/2;
    fNew[n][m][7] = fNew[n][m][5] + (fNew[n][m][2] - fNew[n][m][4])/2;
}

function interiorNodes(n,m){
    fNew[n][m][0] = fOld[n][m][0];
    fNew[n][m][1] = fOld[n][m-1][1];
    fNew[n][m][2] = fOld[n+1][m][2];
    fNew[n][m][3] = fOld[n][m+1][3];
    fNew[n][m][4] = fOld[n-1][m][4];
    fNew[n][m][5] = fOld[n+1][m-1][5];
    fNew[n][m][6] = fOld[n+1][m+1][6];
    fNew[n][m][7] = fOld[n-1][m+1][7];
    fNew[n][m][8] = fOld[n-1][m-1][8];
}

function momentCalculation(n,m){
    rho[n][m] = math.sum(fNew[n][m].slice(0, 9));
    u[n][m] = (fNew[n][m][1] + fNew[n][m][5] + fNew[n][m][8] - fNew[n][m][3] - fNew[n][m][6] - fNew[n][m][7])/rho[n][m];
    v[n][m] = (fNew[n][m][2] + fNew[n][m][5] + fNew[n][m][6] - fNew[n][m][4] - fNew[n][m][7] - fNew[n][m][8])/rho[n][m];
}

function fEqCalculation(n,m,k){
    fEq[n][m][k] = w[k] * rho[n][m] * (1 + dotProduct(Ksi[k], [u[n][m], v[n][m]]) / (cs ** 2) + Math.pow(dotProduct(Ksi[k], [u[n][m], v[n][m]]), 2) / (2 * (cs ** 4)) - (Math.pow(u[n][m], 2) + Math.pow(v[n][m], 2)) / (2 * (cs ** 2)));
}

function collision(){
    for (let i = 0; i < Nx; i++){
        for (let j = 0; j < Ny; j++){
            for (let k = 0; k < 9; k++){
                fOld[i][j][k] = fNew[i][j][k] - (fNew[i][j][k] - fEq[i][j][k]) / Tau;
            }
        }
    }
}

function dotProduct(a, b) {
    return a.reduce((sum, value, index) => sum + value * b[index], 0);
}

// Displaying the Simulation
function setup() {
  pixelDensity(1); // Makes sure the website doesn't try to automatically fill out the area
                   // Without it, the simulation will appear twice on one canvas
  createCanvas(Nx * scaleCell, Ny * scaleCell);
  colorMode(HSB, 255); // Used for hue-based color mapping
  frameRate(60); // Control simulation speed
}

function draw() {
    background(0);

    // for testing resolution
    /*if (frameCount % 60 === 0) {
        console.log(`Frame: ${frameCount}, Time per frame: ${deltaTime.toFixed(2)} ms`);
    }*/

    // Run simulation step
    if (stepCount % drawStep === 0) {
        calculate(); // Data calculation
    }
    stepCount++;
    
    loadPixels(); // create pixels
    // Velocity
    for (let n = 0; n < Nx; n++) {
        for (let m = 0; m < Ny; m++) {
            // For hsb method of color mapping, hue is the direction, red being 0 degrees and blue being 255 degrees
            let velMag = Math.sqrt(u[n][m] ** 2 + v[n][m] ** 2);
            let velAngle = Math.atan2(v[n][m], u[n][m]);
            let hue = map(velAngle, -Math.PI, Math.PI, 0 , 255); // Direction representation
            let brightness = Math.min(255, velMag * scaleColor); // Magnitude representation
             
            let col = color(hue, 255, brightness); 
            let r = red(col);
            let g = green(col);
            let b = blue(col);
            
            for (let dx = 0; dx < scaleCell; dx++) {
                for (let dy = 0; dy < scaleCell; dy++) { 
                    // Display particles    
                    let x = n * scaleCell + dx;
                    let y = (Ny - 1 - m) * scaleCell + dy; // Flip y
                    let index = (x + y * width) * 4;
                    
                    // In p5.js, each pixel is made of 4 elements of the array, i.e. [0] through [3] is the first pixel
                    // Each element represents a RGB value of the pixel
                    pixels[index] = r;       // R
                    pixels[index + 1] = g;   // G
                    pixels[index + 2] = b;   // B
                    pixels[index + 3] = 255; // A (is on a 0 to 255 scale instead of a 0.0 to 1.0 scale)
                    }
                }
        }
    }
    updatePixels(); // update pixels based on coloration after pixels were loaded

    // Pressure
    // Re = (u times d) / [(Tau - 0.5) times cs^2] 
    // Mew = (Tau - 0.5) times cs^2 
}
