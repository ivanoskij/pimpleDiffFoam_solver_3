// Solve the W diffusion equation


//fvScalarMatrix WEqn
//(
//    fvm::ddt(W)
//  + fvm::div(phi, W)
//  - fvm::laplacian(gammaW,W)
// ==
//    fvOptions(W)
//);
//
//WEqn.relax();
//fvOptions.constrain(WEqn);
//WEqn.solve();

////////////////////////////////////////////////


// Solve the W diffusion equation


fvScalarMatrix WEqn
(
    fvm::ddt(W)
  + fvm::div(phi, W)
  + fvm::div(phiW, W)
  - fvm::laplacian(turbulence->nuEff(),W)
 ==
    fvOptions(W)
);

WEqn.relax();
fvOptions.constrain(WEqn);
WEqn.solve();


