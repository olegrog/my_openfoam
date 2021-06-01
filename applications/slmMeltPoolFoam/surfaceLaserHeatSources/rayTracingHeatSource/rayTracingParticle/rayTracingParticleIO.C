/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author(s)
     \\/     M anipulation  |
-------------------------------------------------------------------------------
               DTRMParticle | Copyright (C) 2017-2019 OpenCFD Ltd
         rayTracingParticle | Copyright (C) 2021 Oleg Rogozin
-------------------------------------------------------------------------------
License
    This file is part of slmMeltPoolFoam.

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

#include "rayTracingParticle.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

Foam::string Foam::rayTracingParticle::propertyList_ =
    Foam::rayTracingParticle::propertyList();

const std::size_t Foam::rayTracingParticle::sizeofFields_
(
    sizeof(rayTracingParticle) - sizeof(particle)
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::rayTracingParticle::rayTracingParticle
(
    const polyMesh& mesh,
    Istream& is,
    bool readFields,
    bool newFormat
)
:
    particle(mesh, is, readFields, newFormat),
    p0_(Zero),
    p1_(Zero),
    dQ0_(0),
    dQ_(0),
    isBeingAbsorbed_(-1)
{
    if (readFields)
    {
        if (is.format() == IOstream::ASCII)
        {
            is >> p0_ >> p1_ >> dQ0_ >> dQ_ >> isBeingAbsorbed_;
        }
        else if (!is.checkLabelSize<>() || !is.checkScalarSize<>())
        {
            // Non-native label or scalar size

            is.beginRawRead();

            readRawScalar(is, p0_.data(), vector::nComponents);
            readRawScalar(is, p1_.data(), vector::nComponents);
            readRawScalar(is, &dQ0_);
            readRawScalar(is, &dQ_);
            is.readRaw(reinterpret_cast<char*>(&isBeingAbsorbed_), sizeof(bool));

            is.endRawRead();
        }
        else
        {
            is.read(reinterpret_cast<char*>(&p0_), sizeofFields_);
        }
    }

    is.check(FUNCTION_NAME);
}


void Foam::rayTracingParticle::writeProperties
(
    Ostream& os,
    const wordRes& filters,
    const word& delim,
    const bool namesOnly
) const
{
    particle::writeProperties(os, filters, delim, namesOnly);

    #undef  writeProp
    #define writeProp(Name, Value)                                            \
        particle::writeProperty(os, Name, Value, namesOnly, delim, filters)

    writeProp("p0", p0_);
    writeProp("p1", p1_);
    writeProp("dQ0", dQ0_);
    writeProp("dQ", dQ_);
    writeProp("isBeingAbsorbed", isBeingAbsorbed_);

    #undef writeProp
}


Foam::Ostream& Foam::operator<<(Ostream& os, const rayTracingParticle& p)
{
    if (os.format() == IOstream::ASCII)
    {
        os  << static_cast<const particle&>(p)
            << token::SPACE << p.p0_
            << token::SPACE << p.p1_
            << token::SPACE << p.dQ0_
            << token::SPACE << p.dQ_
            << token::SPACE << p.isBeingAbsorbed_;
    }
    else
    {
        os  << static_cast<const particle&>(p);
        os.write
        (
            reinterpret_cast<const char*>(&p.p0_),
            rayTracingParticle::sizeofFields_
        );
    }

    // Check state of Ostream
    os.check(FUNCTION_NAME);

    return os;
}


// ************************************************************************* //
