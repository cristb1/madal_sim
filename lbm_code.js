// Declaring global variables
const scaleCell = 8; // Size of each cell on screen
const scaleColor = 5; // Scale factor for visualization because velocity magnitudes are really small
const drawStep = 2; // How often to draw frames
let stepCount = 0;

let Nx = 50;
let Ny = 50;
let Re = 100;
let Mew = 0.1;
//let tStop = 1;
    
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

/*function dotProduct(a,b){
    return a.reduce((sum, value, index) => sum + value * b[index], 0);
}*/

// Displaying the Simulation
function setup() {
    pixelDensity(1); // Makes sure the website doesn't try to automatically fill out the area
                   // Without it, the simulation will appear twice on one canvas
    
    // what the sim will be displayed on
    createCanvas(Nx * scaleCell, Ny * scaleCell);
    noSmooth();

    // what the sim itself will be calculated on
    sim = createGraphics(Nx, Ny);
    sim.noSmooth();
    
    frameRate(60); // Control simulation speed
}

function draw() {
    sim.background(0);

    // for testing resolution
    if (frameCount % 60 === 0) {
        console.log(`Frame: ${frameCount}, Time per frame: ${deltaTime.toFixed(2)} ms`);
    }

    // Run simulation step
    if (stepCount % drawStep === 0) {
        // Data calculation
//        for (let t = 0; t < tStop; t++){
        // Streaming for all interior nodes and for boundary nodes
        for (let m = 0; m < Nx; m++){ // x coordinates
            for (let n = 0; n < Ny; n++){ // y coordinates
                if (n == 0){ // top boundary
                    if (m == 0){ // top left corner nodes
                        fNew[n][m][0] = fOld[n][m][0];
                        fNew[n][m][2] = fOld[n+1][m][2];
                        fNew[n][m][3] = fOld[n][m+1][3];
                        fNew[n][m][6] = fOld[n+1][m+1][6];
                        
                        fNew[n][m][1] = fNew[n][m][3];
                        fNew[n][m][4] = fNew[n][m][2];
                        fNew[n][m][8] = fNew[n][m][6];
                        
                        rhoB = (rho[n+1][m] + rho[n][m+1])/2;
                        fNew[n][m][5] = (rhoB - fNew[n][m][0] - fNew[n][m][1] - fNew[n][m][2] - fNew[n][m][3] - fNew[n][m][4] - fNew[n][m][6] - fNew[n][m][8])/2;
                        fNew[n][m][7] = fNew[n][m][5];
                    }else if (m == (Nx - 1)){ // top right corner nodes
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
                }else{ // all other top nodes
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
            }else if(n == (Ny - 1)){ // bottom boundary
                if (m == 0){ // bottom left nodes
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
                }else if (m == (Nx - 1)){ // bottom right nodes
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
                }else{ // all other bottom nodes
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
            }else if (m == 0){ // left boundary
                fNew[n][m][0] = fOld[n][m][0];
                fNew[n][m][2] = fOld[n+1][m][2];
                fNew[n][m][3] = fOld[n][m+1][3];
                fNew[n][m][4] = fOld[n-1][m][4];
                fNew[n][m][6] = fOld[n+1][m+1][6];
                fNew[n][m][7] = fOld[n-1][m+1][7];
                
                fNew[n][m][1] = fNew[n][m][3];
                fNew[n][m][5] = fNew[n][m][7] + (fNew[n][m][4] - fNew[n][m][2])/2;
                fNew[n][m][8] = fNew[n][m][6] + (fNew[n][m][2] - fNew[n][m][4])/2;
            }else if (m == (Nx - 1)){ // right boundary
                fNew[n][m][0] = fOld[n][m][0]; 
                fNew[n][m][1] = fOld[n][m-1][1];
                fNew[n][m][2] = fOld[n+1][m][2];
                fNew[n][m][4] = fOld[n-1][m][4];
                fNew[n][m][5] = fOld[n+1][m-1][5];
                fNew[n][m][8] = fOld[n-1][m-1][8];
                
                fNew[n][m][3] = fNew[n][m][1];
                fNew[n][m][6] = fNew[n][m][8] + (fNew[n][m][4] - fNew[n][m][2])/2;
                fNew[n][m][7] = fNew[n][m][5] + (fNew[n][m][2] - fNew[n][m][4])/2;
            }else{ // interior nodes
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
            
        }

        // moment calculation
        for (let m = 0; m < Nx; m++){
            for (let n = 0; n < Ny; n++){
                let f = fNew[n][m];
                rho[n][m] = f[0] + f[1] + f[2] + f[3] + f[4] + f[5] + f[6] + f[7] + f[8];
                u[n][m] = (fNew[n][m][1] + fNew[n][m][5] + fNew[n][m][8] - fNew[n][m][3] - fNew[n][m][6] - fNew[n][m][7])/rho[n][m];
                v[n][m] = (fNew[n][m][2] + fNew[n][m][5] + fNew[n][m][6] - fNew[n][m][4] - fNew[n][m][7] - fNew[n][m][8])/rho[n][m];
            }
        }

        // fEq calculation
        for (let m = 0; m < Nx; m++){
            for (let n = 0; n < Ny; n++){
                for(let k = 0; k < 9; k++){
                    //fEq[n][m][k] = w[k] * rho[n][m] * (1 + dotProduct(Ksi[k], [u[n][m], v[n][m]]) / (cs ** 2) + Math.pow(dotProduct(Ksi[k], [u[n][m], v[n][m]]), 2) / (2 * (cs ** 4)) - (Math.pow(u[n][m], 2) + Math.pow(v[n][m], 2)) / (2 * (cs ** 2)));
                    let uX = u[n][m];
                    let uY = v[n][m];
                    let KsiX = Ksi[k][0];
                    let KsiY = Ksi[k][1];
                    let KsiDotU = (KsiX * uX) + (KsiY * uY);
                    let uSq = (uX * uX) + (uY * uY);
                    
                    fEq[n][m][k] = w[k] * rho[n][m] * 
                    (1 + 
                    KsiDotU / (cs ** 2) + 
                    (KsiDotU * KsiDotU) / (2 * (cs ** 4)) - 
                    uSq / (2 * (cs ** 2))
                );
                }
            }
        }

        // collision calculation
        for (let i = 0; i < Nx; i++){
            for (let j = 0; j < Ny; j++){
                for (let k = 0; k < 9; k++){
                    fOld[i][j][k] = fNew[i][j][k] - (fNew[i][j][k] - fEq[i][j][k]) / Tau;
                }
            }
        }
    }
//}
    /*}
    stepCount++;*/
//if (stepCount % drawStep === 0) {
    sim.loadPixels(); // create pixels
    // Velocity
    for (let n = 0; n < Nx; n++) {
        for (let m = 0; m < Ny; m++) {
            let velMag = Math.sqrt(u[n][m] ** 2 + v[n][m] ** 2);
            let velMagNormalized = constrain(velMag * scaleColor, 0, 1);
            
            let [r,g,b] = turbo(velMagNormalized);
            // Display particles    
            let x = n;
            let y = (Ny - 1 - m); // Flip y
            let index = (x + y * sim.width) * 4;
                    
            // In p5.js, each pixel is made of 4 elements of the array, i.e. [0] through [3] is the first pixel
            // Each element represents a RGB value of the pixel
            sim.pixels[index] = r;       // R
            sim.pixels[index + 1] = g;   // G
            sim.pixels[index + 2] = b;   // B
            sim.pixels[index + 3] = 255; // A (is on a 0 to 255 scale instead of a 0.0 to 1.0 scale)
            
        }
    }
    sim.updatePixels(); // update pixels based on coloration after pixels were loaded

    background(0);
    image(sim, 0, 0, width, height);
}
stepCount++;

    // Density
    // Re = (u times d) / [(Tau - 0.5) times cs^2] 
    // Mew = (Tau - 0.5) times cs^2 
}
