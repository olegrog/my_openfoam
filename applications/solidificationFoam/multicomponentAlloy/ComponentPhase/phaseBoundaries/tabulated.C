/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author(s)
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2021 Oleg Rogozin
-------------------------------------------------------------------------------
License
    This file is part of solidificationFoam.

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

\*---------------------------------------------------------------------------*/

#include "tabulated.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phaseBoundary::tabulated::tabulated
(
    const dictionary& dict,
    const alloyComponent& component
)
:
    TableBase<scalar>("phaseBoundary", dict)
{
    auto expandedFile(dict.get<fileName>("file"));
    const auto column(dict.get<label>("column"));
    const auto nHeaderLine(dict.getOrDefault<label>("nHeaderLine", 0));
    const char separator(dict.getOrDefault<string>("separator", "\t")[0]);

    autoPtr<ISstream> isPtr(fileHandler().NewIFstream(expandedFile.expand()));
    ISstream& is = *isPtr;

    if (!is.good())
    {
        FatalIOErrorInFunction(is)
            << "Cannot open file " << expandedFile << nl
            << exit(FatalIOError);
    }

    string line;
    label lineNo = 0;
    DynamicList<Tuple2<scalar, scalar>> values;
    DynamicList<string> strings;

    while (is.good())
    {
        is.getLine(line);
        if (lineNo++ < nHeaderLine) continue; // Skip header
        strings.clear();
        std::size_t pos = 0;

        for (label n = 0; (pos != std::string::npos); ++n)
        {
            const auto nPos = line.find(separator, pos);

            if (nPos == std::string::npos)
            {
                strings.append(line.substr(pos));
                pos = nPos;
            }
            else
            {
                strings.append(line.substr(pos, nPos - pos));
                pos = nPos + 1;
            }
        }

        if (strings.size() <= 1)
        {
            break;
        }

        if (strings.size() <= column)
        {
            FatalErrorInFunction
                << "Not enough columns near line " << lineNo
                << ". Require " << (column+1) << " but found "
                << strings << nl
                << exit(FatalError);
        }

        scalar T = readScalar(strings[0]);
        scalar C = readScalar(strings[column]);

        values.append(Tuple2<scalar, scalar>(T, C));
    }

    table_.transfer(values);
}


// ************************************************************************* //
