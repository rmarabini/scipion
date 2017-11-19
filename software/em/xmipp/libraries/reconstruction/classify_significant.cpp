/***************************************************************************
 *
 * Authors:    Carlos Oscar Sanchez Sorzano coss@cnb.csic.es
 *
 * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 * 02111-1307  USA
 *
 *  All comments concerning this program package may be sent to the
 *  e-mail address 'xmipp@cnb.csic.es'
 ***************************************************************************/

#include "classify_significant.h"
#include <data/mask.h>

// Empty constructor =======================================================
ProgClassifySignificant::~ProgClassifySignificant()
{
	for (size_t i=0; i<projector.size(); i++)
		delete projector[i];
}

// Read arguments ==========================================================
void ProgClassifySignificant::readParams()
{
    fnVols = getParam("--ref");
    fnIds = getParam("--id");
    fnAngles = getParam("--angles");
    fnOut = getParam("-o");
    pad = getIntParam("--padding");
}

// Show ====================================================================
void ProgClassifySignificant::show()
{
    if (!verbose)
        return;
    std::cout
    << "Reference volumes:   " << fnVols            << std::endl
    << "Ids:                 " << fnIds             << std::endl
    << "Angles:              " << fnAngles          << std::endl
    << "Output:              " << fnOut             << std::endl
    << "Padding factor:      " << pad               << std::endl
    ;
}

// usage ===================================================================
void ProgClassifySignificant::defineParams()
{
    addUsageLine("Classify a set of images into different classes. See protocol_reconstruct_heterogeneous");
    addParamsLine("   --ref <metadata>            : Reference volumes");
    addParamsLine("   --id <metadata>             : List of itemIds to classified. Sorted.");
    addParamsLine("   --angles <metadata>         : Angles assigned. Each image should have one or several angles");
    addParamsLine("                               : for each volume. The assignment per volume should be organized");
    addParamsLine("                               : in different blocks");
    addParamsLine("   -o <metadata>               : Output metadata with a set of angles per volume");
    addParamsLine("  [--padding <p=2>]            : Padding factor");
}

// Produce side information ================================================
void ProgClassifySignificant::produceSideInfo()
{
	if (verbose>0)
		std::cerr << "Producing side info ..." << std::endl;
    // Read the reference volumes
    Image<double> V;
    MetaData mdVols, mdAngles;
    mdVols.read(fnVols);
    FileName fnVol;
    int i=1;
    FOR_ALL_OBJECTS_IN_METADATA(mdVols)
    {
    	mdVols.getValue(MDL_IMAGE,fnVol,__iter.objId);
    	std::cout << fnVol << std::endl;
        V.read(fnVol);
        V().setXmippOrigin();
        projector.push_back(new FourierProjector(V(),pad,0.5,BSPLINE3));
        currentRowIdx.push_back(0);

        mdAngles.read(formatString("angles_%02d@%s",i,fnAngles.c_str()));
        VMetaData *vmd=new VMetaData();
        mdAngles.asVMetaData(*vmd);
        setAngles.push_back(*vmd);

        subsetAngles.push_back(*(new VMetaData()));
        i+=1;
    }

    // Read the Ids
    MetaData mdIds;
    mdIds.read(fnIds);
    mdIds.getColumnValues(MDL_PARTICLE_ID,setIds);
}

// Current
void ProgClassifySignificant::selectSubset(size_t particleId)
{
	for (size_t i=0; i<projector.size(); i++)
	{
		subsetAngles[i].clear();
		size_t crIdx=currentRowIdx[i];
		MDRow & currentRow=setAngles[i][crIdx];
		size_t currentParticleId;
		currentRow.getValue(MDL_PARTICLE_ID,currentParticleId);
		while (currentParticleId<=particleId)
		{
			subsetAngles[i].push_back(currentRow);
			std::cout << currentRow << std::endl;
			currentRow=setAngles[i][++crIdx];
			currentRow.getValue(MDL_PARTICLE_ID,currentParticleId);
		}
		currentRowIdx[i]=crIdx;
	}
}

//#define DEBUG
void ProgClassifySignificant::run()
{
	show();
	produceSideInfo();

	if (verbose>0)
	{
		std::cerr << "Classifying images ..." << std::endl;
		init_progress_bar(setIds.size());
	}
	for (size_t iid=0; iid<setIds.size(); iid++)
	{
		size_t id=setIds[iid];
		std::cout << id << std::endl;
		selectSubset(id);
		if (verbose>0)
			progress_bar(iid);
	}
	progress_bar(setIds.size());

	/*
    // Read input image and initial parameters
    double rot, tilt, psi;
	rowIn.getValue(MDL_ANGLE_ROT,rot);
	rowIn.getValue(MDL_ANGLE_TILT,tilt);
	rowIn.getValue(MDL_ANGLE_PSI,psi);

	double olda=1.0, oldb=0.0;
	if (rowIn.containsLabel(MDL_CONTINUOUS_GRAY_A)){
		rowIn.getValue(MDL_CONTINUOUS_GRAY_A,olda);
		rowIn.getValue(MDL_CONTINUOUS_GRAY_B,oldb);
	}

	I.read(fnImg);
	I().setXmippOrigin();
	Istddev=I().computeStddev();

    Ifiltered()=I();
    filter.applyMaskSpace(Ifiltered());

	projectVolume(*projector, P, (int)XSIZE(I()), (int)XSIZE(I()),  rot, tilt, psi);
	*/
}
#undef DEBUG
