#   include "readTimeControls.H"

//_AO
//#   include "readPIMPLEControls.H"

    fvSolution solutionDict(runTime);

    const dictionary& pimple = solutionDict.subDict("PIMPLE");

    int nOuterCorr(readInt(pimple.lookup("nOuterCorrectors")));

    int nCorr(readInt(pimple.lookup("nCorrectors")));
    
    int nNonOrthCorr(readInt(pimple.lookup("nNonOrthogonalCorrectors")));

    bool momentumPredictor(readBool(pimple.lookup("momentumPredictor")));    

