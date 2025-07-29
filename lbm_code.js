// This code is what is inside the HTML except that the sliders are integrated differently. In this file, the sliders are dictated through p5.js.
// Since Elementor did not allow that, I had to manually create and style a div for the sliders and its labels when putting the code into the HTML file
// This code is not used directly anymore, but it is the file that I wrote the code on before moving it into the HTML file

// Declaring global variables
const scaleCell = 8; // Size of each cell on screen
const scaleColor = 5; // Scale factor for visualization because velocity magnitudes are really small
const drawStep = 2; // How often to draw frames
let stepCount = 0;
let simTypeSelector;
let velocityFlag = true;
let densityFlag = false;
let mewSlider;
let mewSpan;
let reSlider;
let reSpan;
let gridSlider;
let gridSpan;
let resetFlag = false;

// default values for grid size, Reynolds, and Mew
let Nx = 50;
let Ny = 50;
let Re = 100;
let Mew = 0.1;

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
let L = (Nx - 1)*dx;

let uLid = Re*Mew/L;

let Tau = 3*Mew+0.5;

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

for(let z = 0; z < 9; z++){
for(let m = 0; m < Nx; m++){
    for(let n = 0; n < Ny; n++){
        fOld[n][m][z] = w[z] * fOld[n][m][z];
        fNew[n][m][z] = w[z] * fNew[n][m][z];
        fEq[n][m][z] = w[z] * fEq[n][m][z];
    }
}
}

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
    
    // select type of sim
    simTypeSelector = createSelect();
    simTypeSelector.position(0, Nx * scaleCell + 1);
    simTypeSelector.option('Velocity');
    simTypeSelector.option('Density');
    simTypeSelector.selected('Velocity'); //default visualization type

    // mu slider
    mewSlider = createSlider(0.1,0.32, Mew, 0.01);
    mewSlider.position(150, 400);
    mewSlider.size(80);
    
    mewSpan = createSpan("&#956" + ": " + Mew);
    mewSpan.position(85,404);

    // reynolds number slider
    reSlider = createSlider(1,100, Re);
    reSlider.position(150, 450);
    reSlider.size(80);

    reSpan = createSpan("Re: " + Re);
    reSpan.position(85, 450);

    // Nx slider
    gridSlider = createSlider(50,100, Nx);
    gridSlider.position(150, 500);
    gridSpan = createSpan("Nx: " + Nx);
    gridSpan.position(85, 500);


    noSmooth();

    // what the sim itself will be calculated on
    sim = createGraphics(Nx, Ny);
    sim.noSmooth();
    
    frameRate(60); // Control simulation speed
}

function draw() {
    sim.background(0);
    if(Nx != gridSlider.value() || Mew != (mewSlider.value()) || Re != reSlider.value()){
        resetFlag = true;
    }
    // was here before using the reset flag, will have to make sure removing it doesnt actually break anything before removing
    Mew = (mewSlider.value());
    mewSpan.html("&#956" + ": " + Mew);
    reSpan.html("Re: " + Re);
    gridSpan.html("Nx: " + Nx);
    Re = (reSlider.value());
    Nx = Ny = (gridSlider.value());
    L = (Nx - 1)*dx;
    uLid = Re*Mew/L;
    Tau = 3*Mew+0.5;

// if grid size changes, redo everything   
if(resetFlag){
    // resize canvas
    resizeCanvas(Nx * scaleCell, Ny * scaleCell);
    sim.resizeCanvas(Nx,Ny);

    // resize inputs and text
    simTypeSelector.position(0, Nx * scaleCell + 1);
    mewSlider.position(150, height);
    mewSpan.position(85, height);
    reSlider.position(150, height+50);
    reSpan.position(85, height+50);
    gridSlider.position(150, height+100);
    gridSpan.position(85, height+100);
    
    // reinitialize variables
    Mew = (mewSlider.value());
    mewSpan.html("&#956" + ": " + Mew);
    reSpan.html("Re: " + Re);
    gridSpan.html("Nx: " + Nx);
    Re = (reSlider.value());
    Nx = Ny = (gridSlider.value());
    L = (Nx - 1)*dx;
    uLid = Re*Mew/L;
    Tau = 3*Mew+0.5;    
u = Array.from({length: Nx}, () => 
    Array(Ny).fill(1)
);
v = Array.from({length: Nx}, () => 
    Array(Ny).fill(1)
);
rho = Array.from({length: Nx}, () => 
    Array(Ny).fill(1)
);
fNew = Array.from({length: Nx}, () => 
    Array.from({ length: Ny }, () => 
        Array(9).fill(1)
    )
);
fOld = Array.from({length: Nx}, () => 
    Array.from({ length: Ny }, () => 
        Array(9).fill(1)
    )
);

fEq = Array.from({length: Nx}, () => 
    Array.from({ length: Ny }, () => 
        Array(9).fill(1)
    )
);
for(let z = 0; z < 9; z++){
for(let m = 0; m < Nx; m++){
    for(let n = 0; n < Ny; n++){
        fOld[n][m][z] = w[z] * fOld[n][m][z];
        fNew[n][m][z] = w[z] * fNew[n][m][z];
        fEq[n][m][z] = w[z] * fEq[n][m][z];
    }
}
}
}
resetFlag = false;


    if(simTypeSelector.selected() == 'Velocity'){
        velocityFlag = true;
        densityFlag = false;
    }else if(simTypeSelector.selected() == 'Density'){
        velocityFlag = false;
        densityFlag = true;
    }

    // for testing resolution
    /*if (frameCount % 60 === 0) {
        console.log(`Frame: ${frameCount}, Time per frame: ${deltaTime.toFixed(2)} ms`);
    }*/

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
                        //fNew[n][m][1] = fOld[n][m+1][1]; // added
                        fNew[n][m][2] = fOld[n+1][m][2];
                        fNew[n][m][3] = fOld[n][m+1][3];
                        //fNew[n][m][5] = fOld[n+1][m+1][5]; // added
                        fNew[n][m][6] = fOld[n+1][m+1][6];
                        
                        fNew[n][m][1] = fNew[n][m][3];
                        fNew[n][m][4] = fNew[n][m][2];
                        //fNew[n][m][4] = fNew[n][m][1]; // tried to replace ^^
                        fNew[n][m][8] = fNew[n][m][6];
                        
                        //rhoB = (rho[n+1][m] + rho[n][m+1])/2;
                        rhoB = fNew[n][m].reduce((a,b) => a + b, 0);
                        fNew[n][m][5] = (rhoB - fNew[n][m][0] - fNew[n][m][1] - fNew[n][m][2] - fNew[n][m][3] - fNew[n][m][4] - fNew[n][m][6] - fNew[n][m][8])/2;
                        fNew[n][m][7] = fNew[n][m][5];
                    }else if (m == (Nx - 1)){ // top right corner nodes
                    fNew[n][m][0] = fOld[n][m][0];
                    fNew[n][m][1] = fOld[n][m-1][1];
                    fNew[n][m][2] = fOld[n+1][m][2];
                    //fNew[n][m][3] = fOld[n][m-1][3]; // added
                    fNew[n][m][5] = fOld[n+1][m-1][5];
                    fNew[n][m][6] = fOld[n+1][m][6];
                    
                    fNew[n][m][3] = fNew[n][m][1];
                    fNew[n][m][4] = fNew[n][m][2];
                    //fNew[n][m][4] = fNew[n][m][1]; // tried to replace ^^
                    fNew[n][m][7] = fNew[n][m][5];
                    
                    //rhoB = (rho[n+1][m]+rho[n][m-1])/2;
                    rhoB = fNew[n][m].reduce((a,b) => a + b, 0);
                    fNew[n][m][6] = (rhoB-fNew[n][m][0] - fNew[n][m][1] - fNew[n][m][2] - fNew[n][m][3] - fNew[n][m][4] - fNew[n][m][5] - fNew[n][m][7])/2;
                    //fNew[n][m][5] = (rhoB-fNew[n][m][0] - fNew[n][m][1] - fNew[n][m][2] - fNew[n][m][3] - fNew[n][m][4] - fNew[n][m][6] - fNew[n][m][8])/2;
                    fNew[n][m][8] = fNew[n][m][6];
                }else{ // all other top nodes
                    fNew[n][m][0] = fOld[n][m][0];
                    fNew[n][m][1] = fOld[n][m-1][1];
                    fNew[n][m][2] = fOld[n+1][m][2];
                    fNew[n][m][3] = fOld[n][m+1][3];
                    //fNew[n][m][4] = fOld[n][m+1][4]; // added
                    fNew[n][m][5] = fOld[n+1][m-1][5];
                    fNew[n][m][6] = fOld[n+1][m+1][6];
                    
                    fNew[n][m][4] = fNew[n][m][2];

                    //rhoB = fNew[n][m][0] + fNew[n][m][1] + fNew[n][m][3] + 2*(fNew[n][m][2] + fNew[n][m][5] + fNew[n][m][6]);
                    rhoB = fNew[n][m].reduce((a,b) => a + b, 0);
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
                    //rhoB = (rho[n-1][m] + rho[n][m+1])/2;
                    rhoB = fNew[n][m].reduce((a,b) => a + b, 0);
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
                    //rhoB = (rho[n-1][m] + rho[n][m-1])/2;
                    rhoB = fNew[n][m].reduce((a,b) => a + b, 0);
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
                fNew[n][m][0] = fOld[n][m][0]; // original
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
                //console.log(`f: ${f}`);
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
    if(velocityFlag){
    let minVel = u[0][0];
    let maxVel = u[0][0];
    for(let i = 0; i < Nx; i++){
        for(let j = 0; j < Ny; j++){
            if(maxVel < u[i][j]){
                maxVel = u[i][j];
            }else if(minVel > u[i][j]){
                minVel = u[i][j];
            }
        }
    }
    console.log(`Min is ${minVel}`);    
    for (let m = 0; m < Nx; m++) {
        for (let n = 0; n < Ny; n++) {
            //console.log(u[n][m]);
            let velMag = Math.sqrt(u[n][m] ** 2 + v[n][m] ** 2);
            //let velMagNormalized = constrain(velMag * scaleColor, 0, 1);
            let velMagNormalized = constrain(velMag/uLid, 0, 1);
            
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
}
// Density
if(densityFlag){
    let minDensity = rho[0][0];
    let maxDensity = rho[0][0];
    for(let i = 0; i < Nx; i++){
        for(let j = 0; j < Ny; j++){
            if(maxDensity < rho[i][j]){
                maxDensity = rho[i][j];
            }else if(minDensity > rho[i][j]){
                minDensity = rho[i][j];
            }
        }
    }
    console.log(`Max is ${maxDensity}`);

    for (let n = 0; n < Nx; n++) {
        for (let m = 0; m < Ny; m++) {
            let minDensity = 0.95;
            let maxDensity = 1.1;
            let density = rho[n][m];
            let densityNormalized = constrain(((density - minDensity) / (maxDensity - minDensity)), 0, 1);
            
            let [r,g,b] = turbo(densityNormalized);

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
}
    sim.updatePixels(); // update pixels based on coloration after pixels were loaded

    background(0);
    image(sim, 0, 0, width, height);
}
stepCount++;

}
