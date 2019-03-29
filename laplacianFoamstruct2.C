/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    laplacianFoam

Group
    grpBasicSolvers

Description
    Laplace equation solver for a scalar quantity.

    \heading Solver details
    The solver is applicable to, e.g. for thermal diffusion in a solid.  The
    equation is given by:

    \f[
        \ddt{T}  = \div \left( D_T \grad T \right)
    \f]

    Where:
    \vartable
        T     | Scalar field which is solved for, e.g. temperature
        D_T   | Diffusion coefficient
    \endvartable

    \heading Required fields
    \plaintable
        T     | Scalar field which is solved for, e.g. temperature
    \endplaintable

\*---------------------------------------------------------------------------*/
//#include <iostream>
//#include <fstream>
#include"fvMesh.H"
#include"vector.H"
#include "fvCFD.H"
#include "fvOptions.H"
#include "simpleControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

 //void myDatawrite(labelList abcd,scalarField T,fvMesh& mesh);
 void myDatawrite(labelList abcd,scalarField T,fvMesh& mesh)
{  
 ofstream myFile;
  myFile.open("test.csv",std::ios::app);


    for (int i = 0; i < 4; ++i)
   for (int j=0; j<=2; ++j)
   

    {
        Info<< abcd[i]<< endl;
Info<<mesh.C()[abcd[i]]<<endl;
      myFile<<abcd[i]<<"\t";

      myFile<<T[abcd[i]]<<"\t";
      myFile<<mesh.C()[abcd[i]][j]<<"\n";
    }
 
myFile.close();
}


 

int main(int argc, char *argv[])
{
    #include "addCheckCaseOptions.H"
    #include "setRootCase.H"

    #include "createTime.H"
    #include "createMesh.H"

    simpleControl simple(mesh);

    #include "createFields.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * // 
// write a fun which dumps the value of T on given cells in a text file eg abcd =[1 5 9 13]


    Info<< "\nCalculating temperature distribution\n" << endl;

    while (simple.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        while (simple.correctNonOrthogonal())
        {
            fvScalarMatrix TEqn
            (
                fvm::ddt(T) - fvm::laplacian(DT, T)
             ==
                fvOptions(T)
            );

            fvOptions.constrain(TEqn);
            TEqn.solve();
            fvOptions.correct(T);
        }

        #include "write.H"
 
     labelList abcd(4,0);
abcd[0]=1; 
abcd[1]=5;
abcd[2]=9;
abcd[3]=13;


//fvMesh mesh();
mesh.C()[abcd[0]];
mesh.C()[abcd[1]];
mesh.C()[abcd[2]];
mesh.C()[abcd[3]];

 
     runTime.printExecutionTime(Info);
     myDatawrite(abcd,T,mesh);
      
    }
  
   
    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* /
