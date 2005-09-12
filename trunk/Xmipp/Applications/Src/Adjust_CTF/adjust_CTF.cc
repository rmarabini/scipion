/***************************************************************************
 *
 * Authors:    Carlos Oscar            coss@cnb.uam.es
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
 *  e-mail address 'xmipp@cnb.uam.es'                                  
 ***************************************************************************/

/* INCLUDES ---------------------------------------------------------------- */
#include <Reconstruction/Programs/Prog_adjust_CTF.hh>
#include <XmippData/xmippArgs.hh>

/* MAIN -------------------------------------------------------------------- */
int main (int argc,char *argv[]) {
   Adjust_CTF_Parameters      prog_prm;
   FileName                   fn_in;

// Check the command line
   try {
      fn_in=get_param(argc,argv,"-i");
      prog_prm.read(fn_in);
   } catch (Xmipp_error &XE) {cout << XE; prog_prm.Usage(); exit(1);}

   try {
     // Really adjust CTF
     ROUT_Adjust_CTF(prog_prm);
   } catch (Xmipp_error XE) {cout << XE;}
}

/* Menus ------------------------------------------------------------------- */
