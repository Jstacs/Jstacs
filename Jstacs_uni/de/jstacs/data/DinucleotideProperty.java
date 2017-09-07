/*
 * This file is part of Jstacs.
 * 
 * Jstacs is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later
 * version.
 * 
 * Jstacs is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See the GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License along with
 * Jstacs. If not, see <http://www.gnu.org/licenses/>.
 * 
 * For more information on Jstacs, visit http://www.jstacs.de
 */

package de.jstacs.data;

import java.awt.image.BufferedImage;
import java.util.Arrays;

import de.jstacs.data.alphabets.ContinuousAlphabet;
import de.jstacs.data.alphabets.DNAAlphabet;
import de.jstacs.data.sequences.ArbitrarySequence;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.data.sequences.Sequence.SubSequence;
import de.jstacs.data.sequences.WrongSequenceTypeException;
import de.jstacs.data.sequences.annotation.ReferenceSequenceAnnotation;
import de.jstacs.data.sequences.annotation.SequenceAnnotation;
import de.jstacs.results.ImageResult;
import de.jstacs.utils.REnvironment;

/**
 * This <code>enum</code> defines physicochemical, conformational, and letter-based dinucleotide properties of nucleotide sequences.
 * All dinucleotide parameters are obtained from http://diprodb.fli-leibniz.de/.
 * 
 * @author Jens Keilwagen, Jan Grau
 */
public enum DinucleotideProperty {

	TWIST_KARAS(new double[][]{ { 38.9, 31.12, 32.15, 33.81 },  { 41.41, 34.96, 32.91, 32.15 },  { 41.31, 38.5, 34.96, 31.12 },  { 33.28, 41.31, 41.41, 38.9 } }, "B-DNA", true, "Karas, H. ; Knuppel, R. ; Schulz, W. ; Sklenar, H. ; Wingender, E. (1996): Combining structural analysis of DNA with search routines for the detection of transcription regulatory elements.   Comput. Appl. Biosci. (1996)12 Nr. 5, 441-446", "8996793", HowCreated.CALCULATED, Type.CONFORMATIONAL, "degree", "Values are not listed in the article but referenced to a file that cant be accessed any more. Values taken from PROPERTY DB." ),

	STACKING_ENERGY_SPONER(new double[][]{ { -12, -11.8, -11.5, -10.6 },  { -12.3, -9.5, -13.1, -11.5 },  { -11.4, -13.2, -9.5, -11.8 },  { -11.2, -11.4, -12.3, -12 } }, "B-DNA", true, "Sponer, J. ; Gabb, H. A. ; Leszczynski, J. ; Hobza, P. (1997): Base-Base and Deoxyribose-Base Stacking Interactions in B_DNA and Z-DNA: A Quantum-Chemical Study   Biophysical Journal (1997) 73, 76-87", "9199773", HowCreated.CALCULATED, Type.PHYSICOCHEMICAL, "kcal/mol", "In kcal/mol. Method: Ab initio quantum-chemical method with inclusion of electro correlation. Also given for Z-DNA." ),

	RISE_KARAS(new double[][]{ { 3.16, 3.41, 3.63, 3.89 },  { 3.23, 4.08, 3.6, 3.63 },  { 3.47, 3.81, 4.08, 3.41 },  { 3.21, 3.47, 3.23, 3.16 } }, "B-DNA", true, "Karas, H. ; Knuppel, R. ; Schulz, W. ; Sklenar, H. ; Wingender, E. (1996): Combining structural analysis of DNA with search routines for the detection of transcription regulatory elements.   Comput. Appl. Biosci. (1996)12 Nr. 5, 441-446", "8996793", HowCreated.CALCULATED, Type.CONFORMATIONAL, "angstroem", "Values are not listed in the article but referenced to a file that cant be accessed any more. Values taken from PROPERTY DB.  Calculated by Sklenar, and averaged by Ponomarenko  " ),

	BEND(new double[][]{ { 3.07, 2.97, 2.31, 2.6 },  { 3.58, 2.16, 2.81, 2.31 },  { 2.51, 3.06, 2.16, 2.97 },  { 6.74, 2.51, 3.58, 3.07 } }, "B-DNA", true, "Karas, H. ; Knuppel, R. ; Schulz, W. ; Sklenar, H. ; Wingender, E. (1996): Combining structural analysis of DNA with search routines for the detection of transcription regulatory elements.   Comput. Appl. Biosci. (1996)12 Nr. 5, 441-446", "8996793", HowCreated.CALCULATED, Type.CONFORMATIONAL, "degree", "Values are not listed in the article but referenced to a file that cant be accessed any more. Values taken from PROPERTY DB.  Calculated by Sklenar, and averaged by Ponomarenko" ),

	TIP(new double[][]{ { 1.76, 2, 0.9, 1.87 },  { -1.64, 0.71, 0.22, 0.9 },  { 1.35, 2.5, 0.71, 2 },  { 6.7, 1.35, -1.64, 1.76 } }, "B-DNA", true, "Karas, H. ; Knuppel, R. ; Schulz, W. ; Sklenar, H. ; Wingender, E. (1996): Combining structural analysis of DNA with search routines for the detection of transcription regulatory elements.   Comput. Appl. Biosci. (1996)12 Nr. 5, 441-446", "8996793", HowCreated.CALCULATED, Type.CONFORMATIONAL, "degree", "Values are not listed in the article but referenced to a file that cant be accessed any more. Values taken from PROPERTY DB.  Calculated by Sklenar, and averaged by Ponomarenko" ),

	INCLINATION(new double[][]{ { -1.43, -0.11, -0.92, 0 },  { 1.31, -1.11, 0, 0.92 },  { -0.33, 0, 1.11, 0.11 },  { 0, 0.33, -1.31, 1.43 } }, "DNA", false, "Karas, H. ; Knuppel, R. ; Schulz, W. ; Sklenar, H. ; Wingender, E. (1996): Combining structural analysis of DNA with search routines for the detection of transcription regulatory elements.   Comput. Appl. Biosci. (1996)12 Nr. 5, 441-446", "8996793", HowCreated.CALCULATED, Type.CONFORMATIONAL, "degree", "Values are not listed in the article but referenced to a file that cant be accessed any more. Values taken from PROPERTY DB.  Calculated by Sklenar, and averaged by Ponomarenko" ),

	MAJOR_GROOVE_WIDTH(new double[][]{ { 12.15, 12.37, 13.51, 12.87 },  { 13.58, 15.49, 14.42, 13.51 },  { 13.93, 14.55, 15.49, 12.37 },  { 12.32, 13.93, 13.58, 12.15 } }, "B-DNA", true, "Karas, H. ; Knuppel, R. ; Schulz, W. ; Sklenar, H. ; Wingender, E. (1996): Combining structural analysis of DNA with search routines for the detection of transcription regulatory elements.   Comput. Appl. Biosci. (1996)12 Nr. 5, 441-446", "8996793", HowCreated.CALCULATED, Type.CONFORMATIONAL, "angstroem", "Values are not listed in the article but referenced to a file that cant be accessed any more. Values taken from PROPERTY DB." ),

	MAJOR_GROOVE_DEPTH(new double[][]{ { 9.12, 9.41, 8.96, 8.96 },  { 8.67, 8.45, 8.81, 8.96 },  { 8.76, 8.67, 8.45, 9.41 },  { 9.6, 8.76, 8.67, 9.12 } }, "B-DNA", true, "Karas, H. ; Knuppel, R. ; Schulz, W. ; Sklenar, H. ; Wingender, E. (1996): Combining structural analysis of DNA with search routines for the detection of transcription regulatory elements.   Comput. Appl. Biosci. (1996)12 Nr. 5, 441-446", "8996793", HowCreated.CALCULATED, Type.CONFORMATIONAL, "angstroem", "Values are not listed in the article but referenced to a file that cant be accessed any more. Values taken from PROPERTY DB." ),

	MAJOR_GROOVE_SIZE(new double[][]{ { 3.98, 3.98, 4.7, 4.7 },  { 3.98, 3.98, 4.7, 4.7 },  { 3.26, 3.26, 3.98, 3.98 },  { 3.26, 3.26, 3.98, 3.98 } }, "B-DNA", true, "Gorin, A. A. ; Zhurkin, V. B. ; Olson, W. K. (1995): B-DNA Twisting Correlates with Base-pair Morphology   J. Mol. Biol. (1995) 247, 34-48", "7897660", HowCreated.CALCULATED, Type.CONFORMATIONAL, "angstroem", "" ),

	MAJOR_GROOVE_DISTANCE(new double[][]{ { 3.38, 3.03, 3.36, 3.02 },  { 3.79, 3.38, 3.77, 3.36 },  { 3.4, 3.04, 3.38, 3.03 },  { 3.81, 3.4, 3.79, 3.38 } }, "B-DNA", true, "Gorin, A. A. ; Zhurkin, V. B. ; Olson, W. K. (1995): B-DNA Twisting Correlates with Base-pair Morphology   J. Mol. Biol. (1995) 247, 34-48", "7897660", HowCreated.CALCULATED, Type.CONFORMATIONAL, "angstroem", "" ),

	MINOR_GROOVE_WIDTH(new double[][]{ { 5.3, 6.04, 5.19, 5.31 },  { 4.79, 4.62, 5.16, 5.19 },  { 4.71, 4.74, 4.62, 6.04 },  { 6.4, 4.71, 4.79, 5.3 } }, "B-DNA", true, "Karas, H. ; Knuppel, R. ; Schulz, W. ; Sklenar, H. ; Wingender, E. (1996): Combining structural analysis of DNA with search routines for the detection of transcription regulatory elements.   Comput. Appl. Biosci. (1996)12 Nr. 5, 441-446", "8996793", HowCreated.CALCULATED, Type.CONFORMATIONAL, "angstroem", "Values are not listed in the article but referenced to a file that cant be accessed any more. Values taken from PROPERTY DB." ),

	MINOR_GROOVE_DEPTH(new double[][]{ { 9.03, 8.79, 8.98, 8.91 },  { 9.09, 8.99, 9.06, 8.98 },  { 9.11, 8.98, 8.99, 8.79 },  { 9, 9.11, 9.09, 9.03 } }, "B-DNA", true, "Karas, H. ; Knuppel, R. ; Schulz, W. ; Sklenar, H. ; Wingender, E. (1996): Combining structural analysis of DNA with search routines for the detection of transcription regulatory elements.   Comput. Appl. Biosci. (1996)12 Nr. 5, 441-446", "8996793", HowCreated.CALCULATED, Type.CONFORMATIONAL, "angstroem", "Values are not listed in the article but referenced to a file that cant be accessed any more. Values taken from PROPERTY DB." ),

	MINOR_GROOVE_SIZE(new double[][]{ { 2.98, 3.26, 3.98, 3.26 },  { 3.7, 3.98, 4.7, 3.98 },  { 2.98, 3.26, 3.98, 3.26 },  { 2.7, 2.98, 3.7, 2.98 } }, "B-DNA", true, "Gorin, A. A. ; Zhurkin, V. B. ; Olson, W. K. (1995): B-DNA Twisting Correlates with Base-pair Morphology   J. Mol. Biol. (1995) 247, 34-48", "7897660", HowCreated.CALCULATED, Type.CONFORMATIONAL, "angstroem", "" ),

	MINOR_GROOVE_DISTANCE(new double[][]{ { 2.94, 4.22, 2.79, 4.2 },  { 3.09, 2.8, 3.21, 2.79 },  { 2.95, 4.24, 2.8, 4.22 },  { 2.97, 2.95, 3.09, 2.94 } }, "B-DNA", true, "Gorin, A. A. ; Zhurkin, V. B. ; Olson, W. K. (1995): B-DNA Twisting Correlates with Base-pair Morphology   J. Mol. Biol. (1995) 247, 34-48", "7897660", HowCreated.CALCULATED, Type.CONFORMATIONAL, "angstroem", "" ),

	PERSISTANCE_LENGTH(new double[][]{ { 35, 60, 60, 20 },  { 60, 130, 85, 60 },  { 60, 85, 130, 60 },  { 20, 60, 60, 35 } }, "B-DNA", true, "Hogan, M. E. ; Austin, R. H. (1987): Importance of DNA stiffness in protein-DNA binding specificity   Nature (1987) 329, 263-266", "3627268", HowCreated.EXPERIMENTAL, Type.CONFORMATIONAL, "nanometer", "Describes the rigidity of DNA. Is defined as P = EI / kT. Where E = related to the stress which develops when the long axis of a rod is strained, I = surface moment of inertia for a right cylinder, k = Boltzmann const., T = abs. temp. Values not insi" ),

	MELTING_TEMPERATURE_GOTOH(new double[][]{ { 54.5, 97.73, 58.42, 57.02 },  { 54.71, 85.97, 72.55, 58.42 },  { 86.44, 136.12, 85.97, 97.73 },  { 36.73, 86.44, 54.71, 54.5 } }, "B-DNA", true, "Gotoh, O. ; Tagashira, Y. (1981): Stabilities of Nearest-Neighbor Doublets in Double-Helical DNA Determined by Fitting Calculated Melting Profiles to Observed Profiles   Biopolymers (1981) 20, 1033-1042", "", HowCreated.EXPERIMENTAL, Type.PHYSICOCHEMICAL, "degree", "" ),

	PROBABILITY_CONTACTING_NUCLEOSOME_CORE(new double[][]{ { 18.4, 10.2, 14.5, 7.2 },  { 15.7, 10.2, 1.1, 14.5 },  { 11.3, 5.2, 10.2, 10.2 },  { 6.2, 11.3, 15.7, 18.4 } }, "B-DNA", true, "Hogan, M. E. ; Austin, R. H. (1987): Importance of DNA stiffness in protein-DNA binding specificity   Nature (1987) 329, 263-266", "3627268", HowCreated.CALCULATED, Type.PHYSICOCHEMICAL, "%", "Values are not in the paper. Values taken from the PROPERTY DB." ),

	MOBILITY_TO_BEND_TOWARDS_MAJOR_GROOVE(new double[][]{ { 1.18, 1.06, 1.06, 1.12 },  { 1.06, 0.99, 1.02, 1.04 },  { 1.08, 0.98, 1, 1.02 },  { 1.07, 1.03, 1.03, 1.09 } }, "DNA", false, "Gartenberg, M. R. ; Crothers, D. M. (1988): DNA sequence determinants of CAP-induced bending and protein binding affinity.  Nature (1988)333, 824-829", "2838756", HowCreated.EXPERIMENTAL, Type.PHYSICOCHEMICAL, "mu", "" ),

	MOBILITY_TO_BEND_TOWARDS_MINOR_GROOVE(new double[][]{ { 1.04, 1.1, 1.09, 1.02 },  { 1.16, 1.27, 1.25, 1.16 },  { 1.12, 1.17, 1.25, 1.11 },  { 1.05, 1.2, 1.23, 1.04 } }, "DNA", false, "Gartenberg, M. R. ; Crothers, D. M. (1988): DNA sequence determinants of CAP-induced bending and protein binding affinity.  Nature (1988)333, 824-829", "2838756", HowCreated.EXPERIMENTAL, Type.PHYSICOCHEMICAL, "mu", "" ),

	PROPELLER_TWIST(new double[][]{ { -17.3, -6.7, -14.3, -16.9 },  { -8.6, -12.8, -11.2, -14.3 },  { -15.1, -11.7, -12.8, -6.7 },  { -11.1, -15.1, -8.6, -17.3 } }, "B-DNA", true, "Gorin, A. A. ; Zhurkin, V. B. ; Olson, W. K. (1995): B-DNA Twisting Correlates with Base-pair Morphology   J. Mol. Biol. (1995) 247, 34-48", "7897660", HowCreated.CALCULATED, Type.CONFORMATIONAL, "degree", "" ),

	CLASH_STRENGTH(new double[][]{ { 0.64, 0.95, 2.53, 1.68 },  { 0.8, 1.78, 2.42, 2.53 },  { 0.03, 0.22, 1.78, 0.95 },  { 0, 0.03, 0.8, 0.64 } }, "B-DNA", true, "Gorin, A. A. ; Zhurkin, V. B. ; Olson, W. K. (1995): B-DNA Twisting Correlates with Base-pair Morphology   J. Mol. Biol. (1995) 247, 34-48", "7897660", HowCreated.CALCULATED, Type.CONFORMATIONAL, "angstroem", "Function is defined in the paper, but values are not explicitly listed. Values are taken from the PROPERTY DB." ),

	ENTHALPY_SUGIMOTO(new double[][]{ { -8, -9.4, -6.6, -5.6 },  { -8.2, -10.9, -11.8, -6.6 },  { -8.8, -10.5, -10.9, -9.4 },  { -6.6, -8.8, -8.2, -8 } }, "B-DNA", true, "Sugimoto, N. ; Nakano, S. ; Yoneyama, M. ; Honda, K. (1996): Improved thermodynamic parameters and helix initiation factor to predict stability of DNA duplexes.   Nucleic Acids Research (1996) 24 Nr. 22, 4501-4505", "8948641", HowCreated.EXPERIMENTAL, Type.PHYSICOCHEMICAL, "kcal/mol", "Values given in the paper (delta H)." ),

	ENTROPY_SUGIMOTO(new double[][]{ { -21.9, -25.5, -16.4, -15.2 },  { -21, -28.4, -29, -16.4 },  { -23.5, -26.4, -28.4, -25.5 },  { -18.4, -23.5, -21, -21.9 } }, "B-DNA", true, "Sugimoto, N. ; Nakano, S. ; Yoneyama, M. ; Honda, K. (1996): Improved thermodynamic parameters and helix initiation factor to predict stability of DNA duplexes.   Nucleic Acids Research (1996) 24 Nr. 22, 4501-4505", "8948641", HowCreated.EXPERIMENTAL, Type.PHYSICOCHEMICAL, "cal/mol/K", "Values are given in the paper (delta S)." ),

	SHIFT_RNA(new double[][]{ { -0.08, 0.23, -0.04, -0.06 },  { 0.11, -0.01, 0.3, -0.04 },  { 0.07, 0.07, -0.01, 0.23 },  { -0.02, 0.07, 0.11, -0.08 } }, "A-RNA", true, "Perez, A. ; Noy, A. ; Lankas, F. ; Luque, F. J. ; Orozco, M. (2004): The relative flexibility of B-DNA and A-RNA duplexes: database analysis.   Nucleic Acids Research (2004) 32 No. 20", "15562006", HowCreated.CALCULATED, Type.CONFORMATIONAL, "angstroem", "Values are given in the supplementary." ),

	ROLL_DNA_PROTEIN_COMPLEX_SUZUKI(new double[][]{ { 0.8, -0.2, 5.6, 0 },  { 6.4, 3.3, 6.5, 5.6 },  { 2.4, -2, 3.3, -0.2 },  { 2.7, 2.4, 6.4, 0.8 } }, "B-DNA", true, "Suzuki, M. ; Naoto, Y. ; Finch, J. T. (1996): Role of base-backbone and base-base interactions in alternating DNA conformations   FEBS Letters (1996) 379, 148-152", "8635581", HowCreated.EXPERIMENTAL, Type.CONFORMATIONAL, "degree", "Values are in the Paper under COMPLEX." ),

	TWIST_DNA_PROTEIN_COMPLEX_SUZUKI(new double[][]{ { 35.6, 31.1, 31.9, 29.3 },  { 35.9, 33.3, 34.9, 31.9 },  { 35.9, 34.6, 33.3, 31.1 },  { 39.5, 35.9, 36, 35.6 } }, "B-DNA", true, "Suzuki, M. ; Naoto, Y. ; Finch, J. T. (1996): Role of base-backbone and base-base interactions in alternating DNA conformations   FEBS Letters (1996) 379, 148-152", "8635581", HowCreated.CALCULATED, Type.CONFORMATIONAL, "degree", "Values are in the paper under COMPLEX." ),

	TILT_DNA_PROTEIN_COMPLEX_SUZUKI(new double[][]{ { 1.9, 0.3, 1.3, 0 },  { 0.3, 1, 0, 1.3 },  { 1.7, 0, 1, -0.1 },  { 0, 1.7, 0.3, 1.9 } }, "B-DNA", true, "Suzuki, M. ; Naoto, Y. ; Finch, J. T. (1996): Role of base-backbone and base-base interactions in alternating DNA conformations   FEBS Letters (1996) 379, 148-152", "8635581", HowCreated.CALCULATED, Type.CONFORMATIONAL, "degree", "Values are in the Paper under COMPLEX." ),

	SLIDE_DNA_PROTEIN_COMPLEX_SUZUKI(new double[][]{ { 0.1, -0.6, -0.3, -0.7 },  { 0.4, -0.1, 0.7, -0.3 },  { 0.1, -0.3, -0.1, -0.6 },  { 0.1, 0.1, 0.4, 0.1 } }, "B-DNA", true, "Suzuki, M. ; Naoto, Y. ; Finch, J. T. (1996): Role of base-backbone and base-base interactions in alternating DNA conformations   FEBS Letters (1996) 379, 148-152", "8635581", HowCreated.CALCULATED, Type.CONFORMATIONAL, "angstroem", "Values are in the Paper under COMPLEX." ),

	HYDROPHILICITY_RNA_WEBER(new double[][]{ { 0.023, 0.083, 0.035, 0.09 },  { 0.118, 0.349, 0.193, 0.378 },  { 0.048, 0.146, 0.065, 0.16 },  { 0.112, 0.359, 0.224, 0.389 } }, "RNA", false, "Weber, A. L. ; Lacey, J. C. (1978): Genetic Code Correlations: Amino Acids and Their Anticodon Nucleotides.   J. Mol. Evol. (1978) 11, 199-210", "691071", HowCreated.EXPERIMENTAL, Type.PHYSICOCHEMICAL, "others", "Relative hydrophilicities Rf Values for 16 Dinucleoside Monophosphates in 10/90 v/v 1.0 M ammonium acetate/saturated ammonium sulfate  Original values are in 3'->5' direction. 5'->3' direction can be found in Lacey et al., Org.L.Ev.Bios.(1983)13,3-42" ),

	SHIFT_DNA_PROTEIN_COMPLEX_SUZUKI(new double[][]{ { 0.1, -0.1, -0.2, 0 },  { 0, 0, 0, -0.2 },  { 0.3, 0, 0, -0.1 },  { 0, 0.3, 0, 0.1 } }, "B-DNA", true, "Suzuki, M. ; Naoto, Y. ; Finch, J. T. (1996): Role of base-backbone and base-base interactions in alternating DNA conformations   FEBS Letters (1996) 379, 148-152", "8635581", HowCreated.CALCULATED, Type.CONFORMATIONAL, "angstroem", "Values are in the Paper under COMPLEX." ),

	HYDROPHILICITY_RNA_BARZILAY(new double[][]{ { 0.04, 0.14, 0.08, 0.14 },  { 0.21, 0.49, 0.35, 0.52 },  { 0.1, 0.26, 0.17, 0.27 },  { 0.21, 0.48, 0.34, 0.44 } }, "RNA", false, "Barzilay, I. ; Sussman, J. L. ; Lapidot, Y. (1973): Further Studies on the Chromatographic Behaviour of Dinucleoside Monophosphates   J. Chromatogr. (1973) 79, 139-146", "4350764", HowCreated.EXPERIMENTAL, Type.PHYSICOCHEMICAL, "others", "Rf from paper chromatography (80/18/2:V/V/V).  Original values are in 3'->5' direction. 5'->3' direction can be found in Lacey and Mullins, Orig. Life Evol. Biosph.(1983) 13, 3-42." ),

	RISE_DNA_PROTEIN_COMPLEX_SUZUKI(new double[][]{ { 3.3, 3.4, 3.4, 3.3 },  { 3.4, 3.4, 3.4, 3.4 },  { 3.4, 3.4, 3.4, 3.4 },  { 3.4, 3.4, 3.4, 3.3 } }, "B-DNA", true, "Suzuki, M. ; Naoto, Y. ; Finch, J. T. (1996): Role of base-backbone and base-base interactions in alternating DNA conformations   FEBS Letters (1996) 379, 148-152", "8635581", HowCreated.CALCULATED, Type.CONFORMATIONAL, "angstroem", "Values are in the Paper under COMPLEX." ),

	STACKING_ENERGY_USSERY(new double[][]{ { -5.37, -10.51, -6.78, -6.57 },  { -6.57, -8.26, -9.69, -6.78 },  { -9.81, -14.59, -8.26, -10.51 },  { -3.82, -9.81, -6.57, -5.37 } }, "B-DNA", true, "Ussery, D. W. (2002): DNA Structure: A-, B- and Z-DNA Helix Families   Encyclopedia of Life Sciences  DOI: 10.1038/npg.els.0003122  Article Online Posting Date: May 16, 2002", "", HowCreated.EXPERIMENTAL, Type.PHYSICOCHEMICAL, "kcal/mol", "Values are in the paper but no reference for the original data is given (paper of 1978?)." ),

	FREE_ENERGY_DELCOURT(new double[][]{ { -0.67, -1.28, -1.17, -0.62 },  { -1.19, -1.55, -1.87, -1.17 },  { -1.12, -1.85, -1.55, -1.28 },  { -0.7, -1.12, -1.19, -0.67 } }, "B-DNA", true, "Delcourt, S. G. ; Blake, R. D. (1991): Stacking Energies in DNA   J. of Biol. Chemistry (1991) 266 Nr. 23, 15160-15169", "1869547", HowCreated.EXPERIMENTAL, Type.PHYSICOCHEMICAL, "kcal/mol", "Original and fitted values are given in Frappat et al., Sum rules for free energy and frequency distribution of DNA dinucleotides, Physica A  (2005) 351, 448-460 ." ),

	FREE_ENERGY_BRESLAUER(new double[][]{ { -1.66, -1.13, -1.35, -1.19 },  { -1.8, -2.75, -3.28, -1.35 },  { -1.41, -2.82, -2.75, -1.13 },  { -0.76, -1.41, -1.8, -1.66 } }, "B-DNA", true, "Breslauer, K. J. ; Frank, R. ; Bl?cker, H. ; Marky, L. A. (1986): Predicting DNA duplex stability from the base sequence  Proc. Natl. Acad. Sci. USA (1986) 83, 3746-3750", "3459152", HowCreated.EXPERIMENTAL, Type.PHYSICOCHEMICAL, "kcal/mol", "Original and fitted values are given in Frappat et al., Sum rules for free energy and frequency distribution of DNA dinucleotides, Physica A  (2005) 351, 448-460 ." ),

	FREE_ENERGY_VOLOGODSKII(new double[][]{ { -0.89, -1.35, -1.16, -0.81 },  { -1.37, -1.64, -1.99, -1.16 },  { -1.25, -1.96, -1.64, -1.35 },  { -0.81, -1.16, -1.37, -0.89 } }, "B-DNA", true, "Vologodskii, A. V. ; Amirikyan, B. R. ; Lyubchenko, Y. L. ; Frank-Kamenetskii, M. D. (1984): Allowance for heterogeneous stacking in the DNA helix-coil transition theory.   J. Biomol. Struct. Dyn. (1984) 2, 131-148", "6400927", HowCreated.EXPERIMENTAL, Type.PHYSICOCHEMICAL, "kcal/mol", "Original and fitted values are given in Frappat et al., Sum rules for free energy and frequency distribution of DNA dinucleotides, Physica A  (2005) 351, 448-460 ." ),

	TWIST_DNA_PROTEIN_COMPLEX_OLSON(new double[][]{ { 35.1, 31.5, 31.9, 29.3 },  { 37.3, 32.9, 36.1, 31.9 },  { 36.3, 33.6, 32.9, 31.5 },  { 37.8, 36.3, 37.3, 35.1 } }, "B-DNA", true, "Olson, W. K. ; Gorin, A. A. ; Lu, X. ; Hock, L. M. ; Zhurkin, V. B. (1998): DNA sequence-dependent deformability deduced from protein-DNA crystal complexes   Proc. Natl. Acad. Sci. USA (1998) 95, 11163-11168", "9736707", HowCreated.CALCULATED, Type.CONFORMATIONAL, "degree", "Average values of base pair step parameters computed with CAMPDNA in DNA crystal complexes." ),

	FREE_ENERGY_GOTOH(new double[][]{ { -0.43, -0.98, -0.83, -0.27 },  { -0.97, -1.22, -1.7, -0.83 },  { -0.93, -1.64, -1.22, -0.98 },  { -0.22, -0.93, -0.97, -0.43 } }, "B-DNA", true, "Gotoh, O. ; Tagashira, Y. (1981): Stabilities of Nearest-Neighbor Doublets in Double-Helical DNA Determined by Fitting Calculated Melting Profiles to Observed Profiles   Biopolymers (1981) 20, 1033-1042", "", HowCreated.EXPERIMENTAL, Type.PHYSICOCHEMICAL, "kcal/mol", "Original and fitted values are given in Frappat et al., Sum rules for free energy and frequency distribution of DNA dinucleotides, Physica A  (2005) 351, 448-460 ." ),

	TWIST_TWIST(new double[][]{ { 0.0461, 0.0489, 0.0441, 0.0463 },  { 0.021, 0.0482, 0.0227, 0.0441 },  { 0.0422, 0.0421, 0.0482, 0.0489 },  { 0.0357, 0.0422, 0.021, 0.0461 } }, "B-DNA", true, "Lankas, F. ; Sponer, J. ; Langowski, J. ; Cheatham III, T. E. (2003): DNA Basepair Step Deformability Inferred from Molecular Dynamics Simulations   Biophysical Journal (2003) 85, 2872-2883", "14581192", HowCreated.CALCULATED, Type.CONFORMATIONAL, "kcal/mol degree^2", "Force constants in harmonic potential energy functions describing the deformation of individual basepair steps." ),

	TILT_TILT(new double[][]{ { 0.0389, 0.0411, 0.0371, 0.0404 },  { 0.0275, 0.0414, 0.0278, 0.0371 },  { 0.0392, 0.0396, 0.0414, 0.0411 },  { 0.0245, 0.0392, 0.0275, 0.0389 } }, "B-DNA", true, "Lankas, F. ; Sponer, J. ; Langowski, J. ; Cheatham III, T. E. (2003): DNA Basepair Step Deformability Inferred from Molecular Dynamics Simulations   Biophysical Journal (2003) 85, 2872-2883", "14581192", HowCreated.CALCULATED, Type.CONFORMATIONAL, "kcal/mol degree^2", "Force constants in harmonic potential energy functions describing the deformation of individual basepair steps." ),

	ROLL_ROLL(new double[][]{ { 0.0235, 0.0267, 0.0227, 0.0272 },  { 0.0184, 0.0241, 0.0153, 0.0227 },  { 0.0211, 0.0275, 0.0241, 0.0267 },  { 0.0136, 0.0211, 0.0184, 0.0235 } }, "B-DNA", true, "Lankas, F. ; Sponer, J. ; Langowski, J. ; Cheatham III, T. E. (2003): DNA Basepair Step Deformability Inferred from Molecular Dynamics Simulations   Biophysical Journal (2003) 85, 2872-2883", "14581192", HowCreated.CALCULATED, Type.CONFORMATIONAL, "kcal/mol degree^2", "Force constants in harmonic potential energy functions describing the deformation of individual basepair steps." ),

	TWIST_TILT(new double[][]{ { 0.006, 0.0007, -0.0027, -0.0003 },  { -0.0005, -0.0004, 0.0014, -0.0027 },  { 0.0005, 0.0002, -0.0004, 0.0007 },  { -0.0008, 0.0005, -0.0005, 0.006 } }, "B-DNA", true, "Lankas, F. ; Sponer, J. ; Langowski, J. ; Cheatham III, T. E. (2003): DNA Basepair Step Deformability Inferred from Molecular Dynamics Simulations   Biophysical Journal (2003) 85, 2872-2883", "14581192", HowCreated.CALCULATED, Type.CONFORMATIONAL, "kcal/mol degree^2", "Force constants in harmonic potential energy functions describing the deformation of individual basepair steps." ),

	TWIST_ROLL(new double[][]{ { 0.0083, 0.0076, 0.0057, 0.0081 },  { 0.0049, 0.0044, 0.0031, 0.0057 },  { 0.0086, 0.007, 0.0044, 0.0076 },  { 0.0084, 0.0086, 0.0049, 0.0083 } }, "B-DNA", true, "Lankas, F. ; Sponer, J. ; Langowski, J. ; Cheatham III, T. E. (2003): DNA Basepair Step Deformability Inferred from Molecular Dynamics Simulations   Biophysical Journal (2003) 85, 2872-2883", "14581192", HowCreated.CALCULATED, Type.CONFORMATIONAL, "kcal/mol degree^2", "Force constants in harmonic potential energy functions describing the deformation of individual basepair steps." ),

	TILT_ROLL(new double[][]{ { 0.0033, 0.0029, -0.0027, 0.0007 },  { 0.0009, -0.0009, 0.0011, -0.0027 },  { -0.0002, -0.001, -0.0009, 0.0029 },  { -0.0001, -0.0002, 0.0009, 0.0033 } }, "B-DNA", true, "Lankas, F. ; Sponer, J. ; Langowski, J. ; Cheatham III, T. E. (2003): DNA Basepair Step Deformability Inferred from Molecular Dynamics Simulations   Biophysical Journal (2003) 85, 2872-2883", "14581192", HowCreated.CALCULATED, Type.CONFORMATIONAL, "kcal/mol degree^2", "Force constants in harmonic potential energy functions describing the deformation of individual basepair steps." ),

	SHIFT_SHIFT(new double[][]{ { 1.9748, 1.341, 1.6568, 1.1932 },  { 1.6003, 1.9839, 1.3464, 1.6568 },  { 1.4302, 1.7614, 1.9839, 1.341 },  { 1.5294, 1.4302, 1.6003, 1.9748 } }, "B-DNA", true, "Lankas, F. ; Sponer, J. ; Langowski, J. ; Cheatham III, T. E. (2003): DNA Basepair Step Deformability Inferred from Molecular Dynamics Simulations   Biophysical Journal (2003) 85, 2872-2883", "14581192", HowCreated.CALCULATED, Type.CONFORMATIONAL, "kcal/mol angstroem^2", "Force constants in harmonic potential energy functions describing the deformation of individual basepair steps." ),

	SLIDE_SLIDE(new double[][]{ { 2.9137, 2.9739, 2.7056, 3.3095 },  { 2.2856, 3.2154, 2.0342, 2.7056 },  { 2.5179, 2.7084, 3.2154, 2.9739 },  { 2.2691, 2.5179, 2.2856, 2.9137 } }, "B-DNA", true, "Lankas, F. ; Sponer, J. ; Langowski, J. ; Cheatham III, T. E. (2003): DNA Basepair Step Deformability Inferred from Molecular Dynamics Simulations   Biophysical Journal (2003) 85, 2872-2883", "14581192", HowCreated.CALCULATED, Type.CONFORMATIONAL, "kcal/mol angstroem^2", "Force constants in harmonic potential energy functions describing the deformation of individual basepair steps." ),

	RISE_RISE(new double[][]{ { 7.6206, 9.8821, 6.3875, 10.4992 },  { 6.2903, 7.3347, 4.3896, 6.3875 },  { 8.3295, 10.2808, 7.3347, 9.8821 },  { 5.0546, 8.3295, 6.2903, 7.6206 } }, "B-DNA", true, "Lankas, F. ; Sponer, J. ; Langowski, J. ; Cheatham III, T. E. (2003): DNA Basepair Step Deformability Inferred from Molecular Dynamics Simulations   Biophysical Journal (2003) 85, 2872-2883", "14581192", HowCreated.CALCULATED, Type.CONFORMATIONAL, "kcal/mol angstroem^2", "Force constants in harmonic potential energy functions describing the deformation of individual basepair steps." ),

	SHIFT_SLIDE(new double[][]{ { 0.1711, -0.1574, -0.0263, -0.0965 },  { -0.2832, 0.0572, -0.1867, -0.0263 },  { 0.0259, 0.3178, 0.0572, -0.1574 },  { 0.0516, 0.0259, -0.2832, 0.1711 } }, "B-DNA", true, "Lankas, F. ; Sponer, J. ; Langowski, J. ; Cheatham III, T. E. (2003): DNA Basepair Step Deformability Inferred from Molecular Dynamics Simulations   Biophysical Journal (2003) 85, 2872-2883", "14581192", HowCreated.CALCULATED, Type.CONFORMATIONAL, "kcal/mol anstroem^2", "Force constants in harmonic potential energy functions describing the deformation of individual basepair steps." ),

	SHIFT_RISE(new double[][]{ { 0.1922, -0.0059, -0.0318, -0.0231 },  { -0.0651, 0.2151, -0.0411, -0.0318 },  { 0.025, 0.1312, 0.2151, -0.0059 },  { -0.033, 0.025, -0.0651, 0.1922 } }, "B-DNA", true, "Lankas, F. ; Sponer, J. ; Langowski, J. ; Cheatham III, T. E. (2003): DNA Basepair Step Deformability Inferred from Molecular Dynamics Simulations   Biophysical Journal (2003) 85, 2872-2883", "14581192", HowCreated.CALCULATED, Type.CONFORMATIONAL, "kcal/mol angstroem^2", "Force constants in harmonic potential energy functions describing the deformation of individual basepair steps." ),

	SLIDE_RISE(new double[][]{ { 1.3815, 2.5929, 1.3204, 2.4811 },  { 0.816, 1.1959, 1.4671, 1.3204 },  { 1.1528, 2.5578, 1.1959, 2.5929 },  { 0.913, 1.1528, 0.816, 1.3815 } }, "B-DNA", true, "Lankas, F. ; Sponer, J. ; Langowski, J. ; Cheatham III, T. E. (2003): DNA Basepair Step Deformability Inferred from Molecular Dynamics Simulations   Biophysical Journal (2003) 85, 2872-2883", "14581192", HowCreated.CALCULATED, Type.CONFORMATIONAL, "kcal/mol angstroem^2", "Force constants in harmonic potential energy functions describing the deformation of individual basepair steps." ),

	TWIST_SHIFT(new double[][]{ { 0.0568, 0.0051, -0.0311, -0.0082 },  { -0.0102, 0.0238, 0.0226, -0.0311 },  { -0.0011, -0.0012, 0.0238, 0.0051 },  { -0.0058, -0.0011, -0.0102, 0.0568 } }, "B-DNA", true, "Lankas, F. ; Sponer, J. ; Langowski, J. ; Cheatham III, T. E. (2003): DNA Basepair Step Deformability Inferred from Molecular Dynamics Simulations   Biophysical Journal (2003) 85, 2872-2883", "14581192", HowCreated.CALCULATED, Type.CONFORMATIONAL, "kcal/mol degree angstroem", "Force constants in harmonic potential energy functions describing the deformation of individual basepair steps." ),

	TWIST_SLIDE(new double[][]{ { -0.218, -0.2007, -0.1764, -0.1157 },  { -0.017, -0.225, -0.0855, -0.1764 },  { -0.2056, -0.1929, -0.225, -0.2007 },  { -0.0926, -0.2056, -0.017, -0.218 } }, "B-DNA", true, "Lankas, F. ; Sponer, J. ; Langowski, J. ; Cheatham III, T. E. (2003): DNA Basepair Step Deformability Inferred from Molecular Dynamics Simulations   Biophysical Journal (2003) 85, 2872-2883", "14581192", HowCreated.CALCULATED, Type.CONFORMATIONAL, "kcal/mol degree angstroem", "Force constants in harmonic potential energy functions describing the deformation of individual basepair steps." ),

	TWIST_RISE(new double[][]{ { -0.1587, -0.16, -0.1437, -0.0891 },  { -0.1259, -0.1142, -0.1243, -0.1437 },  { -0.1276, -0.1603, -0.1142, -0.16 },  { -0.0932, -0.1276, -0.1259, -0.1587 } }, "B-DNA", true, "Lankas, F. ; Sponer, J. ; Langowski, J. ; Cheatham III, T. E. (2003): DNA Basepair Step Deformability Inferred from Molecular Dynamics Simulations   Biophysical Journal (2003) 85, 2872-2883", "14581192", HowCreated.CALCULATED, Type.CONFORMATIONAL, "kcal/mol degree angstroem", "Force constants in harmonic potential energy functions describing the deformation of individual basepair steps." ),

	TILT_SHIFT(new double[][]{ { 0.0015, -0.0049, -0.0194, 0.0241 },  { 0.004, -0.0653, -0.0516, -0.0194 },  { -0.0262, -0.0478, -0.0653, -0.0049 },  { 0.0233, -0.0262, 0.004, 0.0015 } }, "B-DNA", true, "Lankas, F. ; Sponer, J. ; Langowski, J. ; Cheatham III, T. E. (2003): DNA Basepair Step Deformability Inferred from Molecular Dynamics Simulations   Biophysical Journal (2003) 85, 2872-2883", "14581192", HowCreated.CALCULATED, Type.CONFORMATIONAL, "kcal/mol degree angstroem", "Force constants in harmonic potential energy functions describing the deformation of individual basepair steps." ),

	TILT_SLIDE(new double[][]{ { -0.0075, -0.0129, 0.0078, -0.0097 },  { -0.0021, 0.005, 0.0103, 0.0078 },  { -0.0023, -0.0183, 0.005, -0.0129 },  { 0.0052, -0.0023, -0.0021, -0.0075 } }, "B-DNA", true, "Lankas, F. ; Sponer, J. ; Langowski, J. ; Cheatham III, T. E. (2003): DNA Basepair Step Deformability Inferred from Molecular Dynamics Simulations   Biophysical Journal (2003) 85, 2872-2883", "14581192", HowCreated.CALCULATED, Type.CONFORMATIONAL, "kcal/mol degree angstroem", "Force constants in harmonic potential energy functions describing the deformation of individual basepair steps." ),

	TILT_RISE(new double[][]{ { -0.2054, 0.0439, 0.0498, 0.0063 },  { -0.0158, -0.0838, 0.0047, 0.0498 },  { -0.0829, -0.0632, -0.0838, 0.0439 },  { -0.0032, -0.0829, -0.0158, -0.2054 } }, "B-DNA", true, "Lankas, F. ; Sponer, J. ; Langowski, J. ; Cheatham III, T. E. (2003): DNA Basepair Step Deformability Inferred from Molecular Dynamics Simulations   Biophysical Journal (2003) 85, 2872-2883", "14581192", HowCreated.CALCULATED, Type.CONFORMATIONAL, "kcal/mol degree angstroem", "Force constants in harmonic potential energy functions describing the deformation of individual basepair steps." ),

	ROLL_SHIFT(new double[][]{ { 0.0158, 0.0141, -0.0143, 0.009 },  { -0.0024, -0.0042, 0.0106, -0.0143 },  { 0.0112, -0.0015, -0.0042, 0.0141 },  { -0.0097, 0.0112, -0.0024, 0.0158 } }, "B-DNA", true, "Lankas, F. ; Sponer, J. ; Langowski, J. ; Cheatham III, T. E. (2003): DNA Basepair Step Deformability Inferred from Molecular Dynamics Simulations   Biophysical Journal (2003) 85, 2872-2883", "14581192", HowCreated.CALCULATED, Type.CONFORMATIONAL, "kcal/mol degree angstroem", "Force constants in harmonic potential energy functions describing the deformation of individual basepair steps." ),

	ROLL_SLIDE(new double[][]{ { -0.022, -0.0022, -0.0291, -0.0499 },  { 0.0093, -0.007, -0.0205, -0.0291 },  { -0.0006, 0.0055, -0.007, -0.0022 },  { -0.0078, -0.0006, 0.0093, -0.022 } }, "B-DNA", true, "Lankas, F. ; Sponer, J. ; Langowski, J. ; Cheatham III, T. E. (2003): DNA Basepair Step Deformability Inferred from Molecular Dynamics Simulations   Biophysical Journal (2003) 85, 2872-2883", "14581192", HowCreated.CALCULATED, Type.CONFORMATIONAL, "kcal/mol degree angstroem", "Force constants in harmonic potential energy functions describing the deformation of individual basepair steps." ),

	ROLL_RISE(new double[][]{ { -0.0541, 0.1089, -0.001, 0.0927 },  { -0.0865, 0.0044, -0.0199, -0.001 },  { -0.0121, 0.1257, 0.0044, 0.1089 },  { -0.037, -0.0121, -0.0865, -0.0541 } }, "B-DNA", true, "Lankas, F. ; Sponer, J. ; Langowski, J. ; Cheatham III, T. E. (2003): DNA Basepair Step Deformability Inferred from Molecular Dynamics Simulations   Biophysical Journal (2003) 85, 2872-2883", "14581192", HowCreated.CALCULATED, Type.CONFORMATIONAL, "kcal/mol degree angstroem", "Force constants in harmonic potential energy functions describing the deformation of individual basepair steps." ),

	STACKING_ENERGY_PEREZ(new double[][]{ { -17.5, -18.1, -15.8, -16.7 },  { -19.5, -14.9, -19.2, -15.8 },  { -14.7, -14.7, -14.9, -18.1 },  { -17, -14.7, -19.5, -17.5 } }, "B-DNA", true, "Perez, A. ; Noy, A. ; Lankas, F. ; Luque, F. J. ; Orozco, M. (2004): The relative flexibility of B-DNA and A-RNA duplexes: database analysis.   Nucleic Acids Research (2004) 32 No. 20", "15562006", HowCreated.CALCULATED, Type.PHYSICOCHEMICAL, "kcal/mol", "Stacking Energy for A-RNA also exists" ),

	TWIST_PEREZ(new double[][]{ { 35, 32, 28, 31 },  { 43, 35, 31, 28 },  { 41, 40, 35, 32 },  { 43, 41, 43, 35 } }, "B-DNA", true, "Perez, A. ; Noy, A. ; Lankas, F. ; Luque, F. J. ; Orozco, M. (2004): The relative flexibility of B-DNA and A-RNA duplexes: database analysis.   Nucleic Acids Research (2004) 32 No. 20", "15562006", HowCreated.CALCULATED, Type.CONFORMATIONAL, "degree", "Values are given in the supplementary." ),

	TILT_PEREZ(new double[][]{ { 0.1, -0.3, 0.2, 0.3 },  { 0, 0.1, 0, 0.2 },  { 0, 0, 0.1, -0.3 },  { -1.4, 0, 0, 0.1 } }, "B-DNA", true, "Perez, A. ; Noy, A. ; Lankas, F. ; Luque, F. J. ; Orozco, M. (2004): The relative flexibility of B-DNA and A-RNA duplexes: database analysis.   Nucleic Acids Research (2004) 32 No. 20", "15562006", HowCreated.CALCULATED, Type.CONFORMATIONAL, "degree", "Values are given in the supplementary." ),

	ROLL_PEREZ(new double[][]{ { 1.4, 1.4, 5.5, -1.2 },  { -1.2, 3.9, 6.2, 5.5 },  { 0.4, -6.8, 3.9, 1.4 },  { -0.6, 0.4, -1.2, 1.4 } }, "B-DNA", true, "Perez, A. ; Noy, A. ; Lankas, F. ; Luque, F. J. ; Orozco, M. (2004): The relative flexibility of B-DNA and A-RNA duplexes: database analysis.   Nucleic Acids Research (2004) 32 No. 20", "15562006", HowCreated.CALCULATED, Type.CONFORMATIONAL, "degree", "Values are given in the supplementary." ),

	SHIFT_PEREZ(new double[][]{ { -0.06, 0.06, 0.06, 0.12 },  { 0.02, 0.05, 0.06, 0.06 },  { 0, -0.3, 0.05, 0.06 },  { -0.17, 0, 0.02, -0.06 } }, "B-DNA", true, "Perez, A. ; Noy, A. ; Lankas, F. ; Luque, F. J. ; Orozco, M. (2004): The relative flexibility of B-DNA and A-RNA duplexes: database analysis.   Nucleic Acids Research (2004) 32 No. 20", "15562006", HowCreated.CALCULATED, Type.CONFORMATIONAL, "angstroem", "Values are given in the supplementary." ),

	SLIDE_PEREZ(new double[][]{ { -0.16, -0.43, 0.34, -0.57 },  { 1.88, 0.28, 0.68, 0.34 },  { -0.01, 0.31, 0.28, -0.43 },  { 0.38, -0.01, 1.88, -0.16 } }, "B-DNA", true, "Perez, A. ; Noy, A. ; Lankas, F. ; Luque, F. J. ; Orozco, M. (2004): The relative flexibility of B-DNA and A-RNA duplexes: database analysis.   Nucleic Acids Research (2004) 32 No. 20", "15562006", HowCreated.CALCULATED, Type.CONFORMATIONAL, "angstroem", "Values are given in the supplementary." ),

	RISE_PEREZ(new double[][]{ { 3.28, 3.23, 3.27, 3.3 },  { 3.32, 3.4, 3.25, 3.27 },  { 3.43, 3.57, 3.4, 3.23 },  { 3.37, 3.43, 3.32, 3.28 } }, "B-DNA", true, "Perez, A. ; Noy, A. ; Lankas, F. ; Luque, F. J. ; Orozco, M. (2004): The relative flexibility of B-DNA and A-RNA duplexes: database analysis.   Nucleic Acids Research (2004) 32 No. 20", "15562006", HowCreated.CALCULATED, Type.CONFORMATIONAL, "angstroem", "Values are given in the supplementary." ),

	SLIDE_STIFFNESS(new double[][]{ { 2.26, 3.03, 2.03, 3.83 },  { 1.78, 1.65, 2, 2.03 },  { 1.93, 2.61, 1.65, 3.03 },  { 1.2, 1.93, 1.78, 2.26 } }, "B-DNA", true, "Goni, J. R. ; Perez, A. ; Torrents, D. ; Orozco, M. (2007): Determining promotor location based on DNA structure first-principles calculations  Genome Biology (2007) 8:R263", "18072969", HowCreated.EXPERIMENTAL, Type.PHYSICOCHEMICAL, "kcal/mol angstroem", "Stiffness constants associated to helical deformations. Used for promoter prediction." ),

	SHIFT_STIFFNESS(new double[][]{ { 1.69, 1.32, 1.46, 1.03 },  { 1.07, 1.43, 1.08, 1.46 },  { 1.32, 1.2, 1.43, 1.32 },  { 0.72, 1.32, 1.07, 1.69 } }, "B-DNA", true, "Goni, J. R. ; Perez, A. ; Torrents, D. ; Orozco, M. (2007): Determining promotor location based on DNA structure first-principles calculations  Genome Biology (2007) 8:R263", "18072969", HowCreated.EXPERIMENTAL, Type.PHYSICOCHEMICAL, "kcal/mol angstroem", "Stiffness constants associated to helical deformations. Used for promoter prediction." ),

	ROLL_STIFFNESS(new double[][]{ { 0.02, 0.023, 0.019, 0.022 },  { 0.017, 0.019, 0.016, 0.019 },  { 0.02, 0.026, 0.019, 0.023 },  { 0.016, 0.02, 0.017, 0.02 } }, "B-DNA", true, "Goni, J. R. ; Perez, A. ; Torrents, D. ; Orozco, M. (2007): Determining promotor location based on DNA structure first-principles calculations  Genome Biology (2007) 8:R263", "18072969", HowCreated.EXPERIMENTAL, Type.PHYSICOCHEMICAL, "kcal/mol degree", "Stiffness constants associated to helical deformations. Used for promoter prediction." ),

	TILT_STIFFNESS(new double[][]{ { 0.038, 0.038, 0.037, 0.036 },  { 0.025, 0.042, 0.026, 0.037 },  { 0.038, 0.036, 0.042, 0.038 },  { 0.018, 0.038, 0.025, 0.038 } }, "B-DNA", true, "Goni, J. R. ; Perez, A. ; Torrents, D. ; Orozco, M. (2007): Determining promotor location based on DNA structure first-principles calculations  Genome Biology (2007) 8:R263", "18072969", HowCreated.EXPERIMENTAL, Type.PHYSICOCHEMICAL, "kcal/mol degree", "Stiffness constants associated to helical deformations. Used for promoter prediction." ),

	TWIST_STIFFNESS(new double[][]{ { 0.026, 0.036, 0.031, 0.033 },  { 0.016, 0.026, 0.014, 0.031 },  { 0.025, 0.025, 0.026, 0.036 },  { 0.017, 0.025, 0.016, 0.026 } }, "B-DNA", true, "Goni, J. R. ; Perez, A. ; Torrents, D. ; Orozco, M. (2007): Determining promotor location based on DNA structure first-principles calculations  Genome Biology (2007) 8:R263", "18072969", HowCreated.EXPERIMENTAL, Type.PHYSICOCHEMICAL, "kcal/mol degree", "Stiffness constants associated to helical deformations. Used for promoter prediction." ),

	FREE_ENERGY_SUGIMOTO(new double[][]{ { -1.2, -1.5, -1.5, -0.9 },  { -1.7, -2.1, -2.8, -1.5 },  { -1.5, -2.3, -2.1, -1.5 },  { -0.9, -1.5, -1.7, -1.2 } }, "B-DNA", true, "Sugimoto, N. ; Nakano, S. ; Yoneyama, M. ; Honda, K. (1996): Improved thermodynamic parameters and helix initiation factor to predict stability of DNA duplexes   Nucleic Acids Research (1996) 24 Nr. 22, 4501-4505", "8948641", HowCreated.EXPERIMENTAL, Type.PHYSICOCHEMICAL, "kcal/mol", "Values are given in the paper (delta G). Fitted values are given in Frappat et al., Sum rules for free energy and frequency distribution of DNA dinucleotides, Physica A  (2005) 351, 448-460 ." ),

	FREE_ENERGY_ALLAWI(new double[][]{ { -1, -1.44, -1.28, -0.88 },  { -1.45, -1.84, -2.17, -1.28 },  { -1.3, -2.24, -1.84, -1.44 },  { -0.58, -1.3, -1.45, -1 } }, "B-DNA", true, "Allawi, H. T. ; SantaLucia, J. (1997): Thermodynamics and NMR of Internal G.T Mismatches in DNA   Biochemistry (1997) 36, 10581-10594", "9265640", HowCreated.EXPERIMENTAL, Type.PHYSICOCHEMICAL, "kcal/mol", "Original and fitted values are given in Frappat et al., Sum rules for free energy and frequency distribution of DNA dinucleotides, Physica A  (2005) 351, 448-460 . (under Unified [18])" ),

	FREE_ENERGY_SANTALUCIA_1996(new double[][]{ { -1.02, -1.43, -1.16, -0.9 },  { -1.7, -1.77, -2.09, -1.16 },  { -1.46, -2.28, -1.77, -1.43 },  { -0.9, -1.46, -1.7, -1.02 } }, "B-DNA", true, "SantaLucia, J. ; Allawi, H. T. ; Seneviratne, P. A. (1996): Improved Nearest-Neighbor Parameters for Predicting DNA Duplex Stability   Biochemistry (1996) 35, 3555-3562", "8639506", HowCreated.EXPERIMENTAL, Type.PHYSICOCHEMICAL, "kcal/mol", "Original and fitted values are given in Frappat et al., Sum rules for free energy and frequency distribution of DNA dinucleotides, Physica A  (2005) 351, 448-460 ." ),

	FREE_ENERGY_OWCZARZY(new double[][]{ { -0.91, -1.25, -1.28, -0.83 },  { -1.54, -1.85, -1.87, -1.28 },  { -1.3, -1.86, -1.85, -1.25 },  { -0.68, -1.3, -1.54, -0.91 } }, "B-DNA", true, "Owczarzy, R. ; Vallone, P. M. ; Goldstein, R. F. ; Benight, A. S. (1992): Studies of DNA Dumbbells VII: Evalualtion of the Next-Nearest-Neighbor Sequence-Dependent Interactions in Duplex DNA   Biopolymers (1992) 32, 29-56", "10737861", HowCreated.EXPERIMENTAL, Type.PHYSICOCHEMICAL, "kcal/mol", "Original and fitted values are given in Frappat et al., Sum rules for free energy and frequency distribution of DNA dinucleotides, Physica A  (2005) 351, 448-460 . (under Benight [6])" ),

	GC_CONTENT(new double[][]{ { 0, 1, 1, 0 },  { 1, 2, 2, 1 },  { 1, 2, 2, 1 },  { 0, 1, 1, 0 } }, "DNA/RNA", false, "Friedel, M. (2008): ", "", HowCreated.CALCULATED, Type.LETTER_BASED, "others", "Each C or G counts +1." ),

	PURINE_AG_CONTENT(new double[][]{ { 2, 1, 2, 1 },  { 1, 0, 1, 0 },  { 2, 1, 2, 1 },  { 1, 0, 1, 0 } }, "DNA/RNA", false, "Friedel, M. (2008): ", "", HowCreated.CALCULATED, Type.LETTER_BASED, "others", "Each purine counts +1." ),

	KETO_GT_CONTENT(new double[][]{ { 0, 0, 0, 1 },  { 0, 0, 1, 1 },  { 1, 1, 2, 2 },  { 1, 1, 1, 2 } }, "DNA/RNA", false, "Friedel, M. (2008): ", "", HowCreated.CALCULATED, Type.LETTER_BASED, "others", "G and T (U) counts +1." ),

	ADENINE_CONTENT(new double[][]{ { 2, 1, 1, 1 },  { 1, 0, 0, 0 },  { 1, 0, 0, 0 },  { 1, 0, 0, 0 } }, "DNA/RNA", false, "Friedel, M. (2008): ", "", HowCreated.CALCULATED, Type.LETTER_BASED, "others", "Each A counts +1." ),

	GUANINE_CONTENT(new double[][]{ { 0, 0, 1, 0 },  { 0, 0, 1, 0 },  { 1, 1, 2, 1 },  { 0, 0, 1, 0 } }, "DNA/RNA", false, "Friedel, M. (2008): ", "", HowCreated.CALCULATED, Type.LETTER_BASED, "others", "Each G counts +1." ),

	CYTOSINE_CONTENT(new double[][]{ { 0, 1, 0, 0 },  { 1, 2, 1, 1 },  { 0, 1, 0, 0 },  { 0, 1, 0, 0 } }, "DNA/RNA", false, "Friedel, M. (2008): ", "", HowCreated.CALCULATED, Type.LETTER_BASED, "others", "Each C counts +1." ),

	THYMINE_CONTENT(new double[][]{ { 0, 0, 0, 1 },  { 0, 0, 0, 1 },  { 0, 0, 0, 1 },  { 1, 1, 1, 2 } }, "DNA/RNA", false, "Friedel, M. (2008): ", "", HowCreated.CALCULATED, Type.LETTER_BASED, "others", "Each T (U) counts +1." ),

	TILT_DNA_PROTEIN_COMPLEX_OLSON(new double[][]{ { -1.4, -0.1, -1.7, 0 },  { 0.5, -0.1, 0, -1.7 },  { -1.5, 0, -0.1, -0.1 },  { 0, -1.5, 0.5, -1.4 } }, "B-DNA", true, "Olson, W. K. ; Gorin, A. A. ; Lu, X. ; Hock, L. M. ; Zhurkin, V. B. (1998): DNA sequence-dependent deformability deduced from protein-DNA crystal complexes   Proc. Natl. Acad. Sci. USA (1998) 95, 11163-11168", "9736707", HowCreated.CALCULATED, Type.CONFORMATIONAL, "angstroem", "Average values of base pair step parameters computed with CAMPDNA in DNA crystal complexes." ),

	ROLL_DNA_PROTEIN_COMPLEX_OLSON(new double[][]{ { 0.7, 0.7, 4.5, 1.1 },  { 4.7, 3.6, 5.4, 4.5 },  { 1.9, 0.3, 3.6, 0.7 },  { 3.3, 1.9, 4.7, 0.7 } }, "B-DNA", true, "Olson, W. K. ; Gorin, A. A. ; Lu, X. ; Hock, L. M. ; Zhurkin, V. B. (1998): DNA sequence-dependent deformability deduced from protein-DNA crystal complexes   Proc. Natl. Acad. Sci. USA (1998) 95, 11163-11168", "9736707", HowCreated.CALCULATED, Type.CONFORMATIONAL, "degree", "Average values of base pair step parameters computed with CAMPDNA in DNA crystal complexes." ),

	SHIFT_DNA_PROTEIN_COMPLEX_OLSON(new double[][]{ { -0.03, 0.13, 0.09, 0 },  { 0.09, 0.05, 0, 0.09 },  { -0.28, 0, 0.05, 0.13 },  { 0, -0.28, 0.09, -0.03 } }, "B-DNA", true, "Olson, W. K. ; Gorin, A. A. ; Lu, X. ; Hock, L. M. ; Zhurkin, V. B. (1998): DNA sequence-dependent deformability deduced from protein-DNA crystal complexes   Proc. Natl. Acad. Sci. USA (1998) 95, 11163-11168", "9736707", HowCreated.CALCULATED, Type.CONFORMATIONAL, "angstroem", "Average values of base pair step parameters computed with CAMPDNA in DNA crystal complexes." ),

	SLIDE_DNA_PROTEIN_COMPLEX_OLSON(new double[][]{ { -0.08, -0.58, -0.25, -0.59 },  { 0.53, -0.22, 0.41, -0.25 },  { 0.09, -0.38, -0.22, -0.58 },  { 0.05, 0.09, 0.53, -0.08 } }, "B-DNA", true, "Olson, W. K. ; Gorin, A. A. ; Lu, X. ; Hock, L. M. ; Zhurkin, V. B. (1998): DNA sequence-dependent deformability deduced from protein-DNA crystal complexes   Proc. Natl. Acad. Sci. USA (1998) 95, 11163-11168", "9736707", HowCreated.CALCULATED, Type.CONFORMATIONAL, "angstroem", "Average values of base pair step parameters computed with CAMPDNA in DNA crystal complexes." ),

	RISE_DNA_PROTEIN_COMPLEX_OLSON(new double[][]{ { 3.27, 3.36, 3.34, 3.31 },  { 3.33, 3.42, 3.39, 3.34 },  { 3.37, 3.4, 3.42, 3.36 },  { 3.42, 3.37, 3.33, 3.27 } }, "B-DNA", true, "Olson, W. K. ; Gorin, A. A. ; Lu, X. ; Hock, L. M. ; Zhurkin, V. B. (1998): DNA sequence-dependent deformability deduced from protein-DNA crystal complexes   Proc. Natl. Acad. Sci. USA (1998) 95, 11163-11168", "9736707", HowCreated.CALCULATED, Type.CONFORMATIONAL, "angstroem", "Average values of base pair step parameters computed with CAMPDNA in DNA crystal complexes." ),

	TWIST_GORIN(new double[][]{ { 35.8, 35.8, 30.5, 33.4 },  { 36.9, 33.4, 31.1, 30.5 },  { 39.3, 38.3, 33.4, 35.8 },  { 40, 39.3, 36.9, 35.8 } }, "B-DNA", true, "Gorin, A. A. ; Zhurkin, V. B. ; Olson, W. K. (1995): B-DNA Twisting Correlates with Base-pair Morphology   J. Mol. Biol. (1995) 247, 34-48", "7897660", HowCreated.EXPERIMENTAL, Type.CONFORMATIONAL, "degree", "Rotation around the helix axis (Wikipedia)." ),

	TILT_GORIN(new double[][]{ { -0.4, -0.9, -2.6, 0 },  { 0.6, -1.1, 0, -2.6 },  { -0.4, 0, -1.1, -0.9 },  { 0, -0.4, 0.6, -0.4 } }, "B-DNA", true, "Gorin, A. A. ; Zhurkin, V. B. ; Olson, W. K. (1995): B-DNA Twisting Correlates with Base-pair Morphology   J. Mol. Biol. (1995) 247, 34-48", "7897660", HowCreated.CALCULATED, Type.CONFORMATIONAL, "degree", "Rotation around an axis in the base-pair plane perpendicular to the first." ),

	ROLL_GORIN(new double[][]{ { 0.5, 0.4, 2.9, -0.6 },  { 1.1, 6.5, 6.6, 2.9 },  { -0.1, -7, 6.5, 0.4 },  { 2.6, -0.1, 1.1, 0.5 } }, "B-DNA", true, "Gorin, A. A. ; Zhurkin, V. B. ; Olson, W. K. (1995): B-DNA Twisting Correlates with Base-pair Morphology   J. Mol. Biol. (1995) 247, 34-48", "7897660", HowCreated.CALCULATED, Type.CONFORMATIONAL, "degree", "Rotation around an axis in the plane of the base pair directed from one strand to the other (Wikipedia)." ),

	SLIDE_GORIN(new double[][]{ { -0.03, -0.13, 0.47, -0.37 },  { 1.46, 0.6, 0.63, 0.47 },  { -0.07, 0.29, 0.6, -0.13 },  { 0.74, -0.07, 1.46, -0.03 } }, "B-DNA", true, "Gorin, A. A. ; Zhurkin, V. B. ; Olson, W. K. (1995): B-DNA Twisting Correlates with Base-pair Morphology   J. Mol. Biol. (1995) 247, 34-48", "7897660", HowCreated.CALCULATED, Type.CONFORMATIONAL, "angstroem", "Displacement along an axis in the plane of the base pair directed from one strand to the other." ),

	TWIST_SUZUKI(new double[][]{ { 35.3, 32.6, 31.2, 31.2 },  { 39.2, 33.3, 36.6, 31.2 },  { 40.3, 37.3, 33.3, 32.6 },  { 40.5, 40.3, 39.2, 35.3 } }, "B-DNA", true, "Suzuki, M. ; Naoto, Y. ; Finch, J. T. (1996): Role of base-backbone and base-base interactions in alternating DNA conformations   FEBS Letters (1996) 379, 148-152", "8635581", HowCreated.CALCULATED, Type.CONFORMATIONAL, "degree", "Rotation around the helix axis (Wikipedia)." ),

	TILT_SUZUKI(new double[][]{ { 0.5, 0.1, 2.8, 0 },  { -0.7, 2.7, 0, 2.8 },  { 0.9, 0, 2.7, 0.1 },  { 0, 0.9, -0.7, 0.5 } }, "B-DNA", true, "Suzuki, M. ; Naoto, Y. ; Finch, J. T. (1996): Role of base-backbone and base-base interactions in alternating DNA conformations   FEBS Letters (1996) 379, 148-152", "8635581", HowCreated.CALCULATED, Type.CONFORMATIONAL, "degree", "Rotation around an axis in the base-pair plane perpendicular to the first." ),

	ROLL_SUZUKI(new double[][]{ { 0.3, 0.5, 4.5, -0.8 },  { 0.5, 6, 3.1, 4.5 },  { -1.3, -6.2, 6, 0.5 },  { 2.8, -1.3, 0.5, 0.3 } }, "B-DNA", true, "Suzuki, M. ; Naoto, Y. ; Finch, J. T. (1996): Role of base-backbone and base-base interactions in alternating DNA conformations   FEBS Letters (1996) 379, 148-152", "8635581", HowCreated.CALCULATED, Type.CONFORMATIONAL, "degree", "Rotation around an axis in the plane of the base pair directed from one strand to the other (Wikipedia)." ),

	SHIFT_SUZUKI(new double[][]{ { 0, 0.2, -0.4, 0 },  { 0.1, 0, 0, -0.4 },  { 0, 0, 0, 0.2 },  { 0, 0, 0.1, 0 } }, "B-DNA", true, "Suzuki, M. ; Naoto, Y. ; Finch, J. T. (1996): Role of base-backbone and base-base interactions in alternating DNA conformations   FEBS Letters (1996) 379, 148-152", "8635581", HowCreated.CALCULATED, Type.CONFORMATIONAL, "angstroem", "Def.: Displacement along an axis in the base-pair plane perpendicular to the first, directed from the minor to the major groove." ),

	SLIDE_SUZUKI(new double[][]{ { -0.1, -0.2, 0.4, -0.4 },  { 1.6, 0.8, 0.7, 0.4 },  { 0, 0.4, 0.8, -0.2 },  { 0.9, 0, 1.6, -0.1 } }, "B-DNA", true, "Suzuki, M. ; Naoto, Y. ; Finch, J. T. (1996): Role of base-backbone and base-base interactions in alternating DNA conformations   FEBS Letters (1996) 379, 148-152", "8635581", HowCreated.CALCULATED, Type.CONFORMATIONAL, "angstroem", "Displacement along an axis in the plane of the base pair directed from one strand to other." ),

	RISE_SUZUKI(new double[][]{ { 3.3, 3.3, 3.3, 3.3 },  { 3.4, 3.4, 3.4, 3.3 },  { 3.3, 3.5, 3.4, 3.3 },  { 3.4, 3.3, 3.4, 3.3 } }, "B-DNA", true, "Suzuki, M. ; Naoto, Y. ; Finch, J. T. (1996): Role of base-backbone and base-base interactions in alternating DNA conformations   FEBS Letters (1996) 379, 148-152", "8635581", HowCreated.CALCULATED, Type.CONFORMATIONAL, "angstroem", "Displacement along the helix axis (Wikipedia)." ),

	TWIST_SHPIGELMAN(new double[][]{ { 35.62, 34.4, 27.7, 31.5 },  { 34.5, 33.67, 29.8, 27.7 },  { 36.9, 40, 33.67, 34.4 },  { 36, 36.9, 34.5, 35.62 } }, "B-DNA", true, "Shpigelman, E. S. ; Trifonov, E. N. ; Bolshoy, A. (1993): CURVATURE: software for the analysis of curved DNA   Bioinformatics (1993) 9 Nr. 4, 435-440", "8402210", HowCreated.CALCULATED, Type.CONFORMATIONAL, "degree", "Rotation around the helix axis (Wikipedia).  Data taken from Kabsch et al.: NAR (1982) 10, 1097-1104" ),

	WEDGE(new double[][]{ { 7.2, 1.1, 8.4, 2.6 },  { 3.5, 2.1, 6.7, 8.4 },  { 5.3, 5, 2.1, 1.1 },  { 0.9, 5.3, 3.5, 7.2 } }, "B-DNA", true, "Shpigelman, E. S. ; Trifonov, E. N. ; Bolshoy, A. (1993): CURVATURE: software for the analysis of curved DNA   Bioinformatics (1993) 9 Nr. 4, 435-440", "8402210", HowCreated.CALCULATED, Type.CONFORMATIONAL, "degree", "Data taken from Bolshoy et al.: PNAS (1991) 88, 2312-2316" ),

	DIRECTION(new double[][]{ { -154, 143, 2, 0 },  { -64, -57, 0, -2 },  { 120, 180, 57, -143 },  { 0, -120, 64, 154 } }, "DNA", false, "Shpigelman, E. S. ; Trifonov, E. N. ; Bolshoy, A. (1993): CURVATURE: software for the analysis of curved DNA   Bioinformatics (1993) 9 Nr. 4, 435-440", "8402210", HowCreated.CALCULATED, Type.CONFORMATIONAL, "degree", "Data taken from Bolshoy et al.: PNAS (1991) 88, 2312-2316  Describes the direction of the deflection angle (wedge angle)." ),

	SLIDE_RNA(new double[][]{ { -1.27, -1.43, -1.5, -1.36 },  { -1.46, -1.78, -1.89, -1.5 },  { -1.7, -1.39, -1.78, -1.43 },  { -1.45, -1.7, -1.46, -1.27 } }, "A-RNA", true, "Perez, A. ; Noy, A. ; Lankas, F. ; Luque, F. J. ; Orozco, M. (2004): The relative flexibility of B-DNA and A-RNA duplexes: database analysis.   Nucleic Acids Research (2004) 32 No. 20", "15562006", HowCreated.CALCULATED, Type.CONFORMATIONAL, "angstroem", "Values are given in the supplementary." ),

	RISE_RNA(new double[][]{ { 3.18, 3.24, 3.3, 3.24 },  { 3.09, 3.32, 3.3, 3.3 },  { 3.38, 3.22, 3.32, 3.24 },  { 3.26, 3.38, 3.09, 3.18 } }, "A-RNA", true, "Perez, A. ; Noy, A. ; Lankas, F. ; Luque, F. J. ; Orozco, M. (2004): The relative flexibility of B-DNA and A-RNA duplexes: database analysis.   Nucleic Acids Research (2004) 32 No. 20", "15562006", HowCreated.CALCULATED, Type.CONFORMATIONAL, "angstroem", "Values are given in the supplementary." ),

	TILT_RNA(new double[][]{ { -0.8, 0.8, 0.5, 1.1 },  { 1, 0.3, -0.1, 0.5 },  { 1.3, 0, 0.3, 0.8 },  { -0.2, 1.3, 1, -0.8 } }, "A-RNA", true, "Perez, A. ; Noy, A. ; Lankas, F. ; Luque, F. J. ; Orozco, M. (2004): The relative flexibility of B-DNA and A-RNA duplexes: database analysis.   Nucleic Acids Research (2004) 32 No. 20", "15562006", HowCreated.CALCULATED, Type.CONFORMATIONAL, "degree", "Values are given in the supplementary." ),

	ROLL_RNA(new double[][]{ { 7, 4.8, 8.5, 7.1 },  { 9.9, 8.7, 12.1, 8.5 },  { 9.4, 6.1, 12.1, 4.8 },  { 10.7, 9.4, 9.9, 7 } }, "A-RNA", true, "Perez, A. ; Noy, A. ; Lankas, F. ; Luque, F. J. ; Orozco, M. (2004): The relative flexibility of B-DNA and A-RNA duplexes: database analysis.   Nucleic Acids Research (2004) 32 No. 20", "15562006", HowCreated.CALCULATED, Type.CONFORMATIONAL, "degree", "Values are given in the supplementary." ),

	TWIST_RNA(new double[][]{ { 31, 32, 30, 33 },  { 31, 32, 27, 30 },  { 32, 35, 32, 32 },  { 32, 32, 31, 31 } }, "A-RNA", true, "Perez, A. ; Noy, A. ; Lankas, F. ; Luque, F. J. ; Orozco, M. (2004): The relative flexibility of B-DNA and A-RNA duplexes: database analysis.   Nucleic Acids Research (2004) 32 No. 20", "15562006", HowCreated.CALCULATED, Type.CONFORMATIONAL, "degree", "Values are given in the supplementary." ),

	STACKING_ENERGY_RNA(new double[][]{ { -13.7, -13.8, -14, -15.4 },  { -14.4, -11.1, -15.6, -14 },  { -14.2, -16.9, -11.1, -13.8 },  { -16, -14.2, -14.4, -13.7 } }, "A-RNA", true, "Perez, A. ; Noy, A. ; Lankas, F. ; Luque, F. J. ; Orozco, M. (2004): The relative flexibility of B-DNA and A-RNA duplexes: database analysis.   Nucleic Acids Research (2004) 32 No. 20", "15562006", HowCreated.CALCULATED, Type.PHYSICOCHEMICAL, "kcal/mol", "Values given in the paper." ),

	RISE_STIFFNESS(new double[][]{ { 7.65, 8.93, 7.08, 9.07 },  { 6.38, 8.04, 6.23, 7.08 },  { 8.56, 9.53, 8.04, 8.93 },  { 6.23, 8.56, 6.38, 7.65 } }, "B-DNA", true, "Goni, J. R. ; Perez, A. ; Torrents, D. ; Orozco, M. (2007): Determining promotor location based on DNA structure first-principles calculations  Genome Biology (2007) 8:R263", "18072969", HowCreated.EXPERIMENTAL, Type.PHYSICOCHEMICAL, "kcal/mol angstroem", "Stiffness constants associated to helical deformations. Used for promoter prediction." ),

	MELTING_TEMPERATURE_ANSELMI(new double[][]{ { 0.945, 1.07, 0.956, 0.952 },  { 0.945, 1.036, 0.997, 0.956 },  { 1.037, 1.18, 1.036, 1.07 },  { 0.894, 1.037, 0.945, 0.945 } }, "B-DNA", true, "Anselmi, C. ; Santis, P. D. ; Paparcone, R. ; Savino, M. ; Scipioni, A. (2002): From the sequence to the superstructural properties of DNAs.  Biophysical Chemistry (2002) 95, 23-47", "11880171", HowCreated.EXPERIMENTAL, Type.PHYSICOCHEMICAL, "others", "Normalized melting temperatures of O. Gotoh and Y. Tagashira (1981), Biopolymers, 20, 1033-1042" ),

	STACKING_ENERGY_ANSELMI(new double[][]{ { 0.703, 1.323, 0.78, 0.854 },  { 0.79, 0.984, 1.124, 0.78 },  { 1.23, 1.792, 0.984, 1.323 },  { 0.615, 1.23, 0.79, 0.703 } }, "B-DNA", true, "Anselmi, C. ; Santis, P. D. ; Paparcone, R. ; Savino, M. ; Scipioni, A. (2002): From the sequence to the superstructural properties of DNAs.  Biophysical Chemistry (2002) 95, 23-47", "11880171", HowCreated.EXPERIMENTAL, Type.PHYSICOCHEMICAL, "others", "normalized values from R.L. Ornsteinet al., An  optimized potential function for the calculation of  nucleic acid interaction energies. I. Base stacking, Biopolymers 17 (1978) 2341?2360." ),

	ENTHALPY_RNA_FREIER(new double[][]{ { -6.6, -10.2, -7.6, -5.7 },  { -10.5, -12.2, -8, -7.6 },  { -13.3, -14.2, -12.2, -10.2 },  { -8.1, -10.2, -7.6, -6.6 } }, "A-RNA", true, "Freier, S. M. ; Kierzek, R. ; Jaeger, J. A. ; Sugimoto, N. ; Caruthers, M. H. ; Neilson, T. ; Turner, D. H. (1986): Improved free-energy parameters for predictions of RNA duplex stability  Proc. Natl. Acad. Sci. USA (1986) 83, 9373-9377", "2432595", HowCreated.EXPERIMENTAL, Type.PHYSICOCHEMICAL, "kcal/mol", "In the paper under delta H." ),

	ENTROPY_RNA_FREIER(new double[][]{ { -18.4, -26.2, -19.2, -15.5 },  { -27.8, -29.7, -19.4, -19.2 },  { -35.5, -34.9, -29.7, -26.2 },  { -22.6, -26.2, -19.2, -18.4 } }, "A-RNA", true, "Freier, S. M. ; Kierzek, R. ; Jaeger, J. A. ; Sugimoto, N. ; Caruthers, M. H. ; Neilson, T. ; Turner, D. H. (1986): Improved free-energy parameters for predictions of RNA duplex stability  Proc. Natl. Acad. Sci. USA (1986) 83, 9373-9377", "2432595", HowCreated.EXPERIMENTAL, Type.PHYSICOCHEMICAL, "eu", "In the Paper under delta S [eu]." ),

	FREE_ENERGY_RNA_FREIER(new double[][]{ { -0.9, -2.1, -1.7, -0.9 },  { -1.8, -2.9, -2, -1.7 },  { -2.3, -3.4, -2.9, -2.1 },  { -1.1, -2.1, -1.7, -0.9 } }, "A-RNA", true, "Freier, S. M. ; Kierzek, R. ; Jaeger, J. A. ; Sugimoto, N. ; Caruthers, M. H. ; Neilson, T. ; Turner, D. H. (1986): Improved free-energy parameters for predictions of RNA duplex stability  Proc. Natl. Acad. Sci. USA (1986) 83, 9373-9377", "2432595", HowCreated.EXPERIMENTAL, Type.PHYSICOCHEMICAL, "kcal/mol", "In the paper under delta G." ),

	FREE_ENERGY_RNA_XIA(new double[][]{ { -0.93, -2.24, -2.08, -1.1 },  { -2.11, -3.26, -2.36, -2.08 },  { -2.35, -3.42, -3.26, -2.24 },  { -1.33, -2.35, -2.11, -0.93 } }, "A-RNA", true, "Xia, T.; SantaLucia, J.; Burkard, M. E.; Kierzek, R.; Schroeder, S. J.; Jiao, X; Cox, C.; Turner, D. H. (1998): Thermodynamic Parameters for an Expanded Nearest-Neighbor TrainableStatisticalModel for Formation of RNA Duplexes with Watson-Crick Base Pairs  Biochemistry (1998) 37, 14719-14735", "9778347", HowCreated.EXPERIMENTAL, Type.PHYSICOCHEMICAL, "kcal/mol", "Table 4, delta G" ),

	ENTHALPY_RNA_XIA(new double[][]{ { -6.82, -11.4, -10.48, -9.38 },  { -10.44, -13.39, -10.64, -10.48 },  { -12.44, -14.88, -13.39, -11.4 },  { -7.69, -12.44, -10.44, -6.82 } }, "A-RNA", true, "Xia, T.; SantaLucia, J.; Burkard, M. E.; Kierzek, R.; Schroeder, S. J.; Jiao, X; Cox, C.; Turner, D. H. (1998): Thermodynamic Parameters for an Expanded Nearest-Neighbor TrainableStatisticalModel for Formation of RNA Duplexes with Watson-Crick Base Pairs  Biochemistry (1998) 37, 14719-14735", "9778347", HowCreated.EXPERIMENTAL, Type.PHYSICOCHEMICAL, "kcal/mol", "Table 4, delta H" ),

	ENTROPY_RNA_XIA(new double[][]{ { -19, -29.5, -27.1, -26.7 },  { -26.9, -32.7, -26.7, -27.1 },  { -32.5, -36.9, -32.7, -29.5 },  { -20.5, -32.5, -26.9, -19 } }, "A-RNA", true, "Xia, T.; SantaLucia, J.; Burkard, M. E.; Kierzek, R.; Schroeder, S. J.; Jiao, X; Cox, C.; Turner, D. H. (1998): Thermodynamic Parameters for an Expanded Nearest-Neighbor TrainableStatisticalModel for Formation of RNA Duplexes with Watson-Crick Base Pairs  Biochemistry (1998) 37, 14719-14735", "9778347", HowCreated.EXPERIMENTAL, Type.PHYSICOCHEMICAL, "eu", "[eu] Table 4, delta G" ),

	ROLL_ANSELMI(new double[][]{ { -5.4, -2.5, 1, -7.3 },  { 6.8, 1.3, 4.6, 1 },  { 2, -3.7, 1.3, -2.5 },  { 8, 2, 6.8, -5.4 } }, "B-DNA", false, "Anselmi, C. ; Santis, P. D. ; Paparcone, R. ; Savino, M. ; Scipioni, A. (2002): From the sequence to the superstructural properties of DNAs.  Biophysical Chemistry (2002) 95, 23-47", "11880171", HowCreated.CALCULATED, Type.CONFORMATIONAL, "Dimension", "" ),

	TILT_ANSELMI(new double[][]{ { -0.5, -2.7, -1.6, 0 },  { 0.4, 0.6, 0, 1.6 },  { -1.7, 0, -0.6, 2.7 },  { 0, 1.7, -0.4, 0.5 } }, "B-DNA", false, "Anselmi, C. ; Santis, P. D. ; Paparcone, R. ; Savino, M. ; Scipioni, A. (2002): From the sequence to the superstructural properties of DNAs.  Biophysical Chemistry (2002) 95, 23-47", "11880171", HowCreated.CALCULATED, Type.CONFORMATIONAL, "Dimension", "" ),

	TWIST_ANSELMI(new double[][]{ { 36, 33.7, 34.4, 35.3 },  { 34.1, 33.1, 33.5, 34.4 },  { 34.6, 33.3, 33.1, 33.7 },  { 34.5, 34.6, 34.1, 36 } }, "B-DNA", false, "Anselmi, C. ; Santis, P. D. ; Paparcone, R. ; Savino, M. ; Scipioni, A. (2002): From the sequence to the superstructural properties of DNAs.  Biophysical Chemistry (2002) 95, 23-47", "11880171", HowCreated.CALCULATED, Type.CONFORMATIONAL, "Dimension", "" ),

	ROLL_PACKER(new double[][]{ { 2.3, -2, 0.5, -8.1 },  { 7.4, 1.4, 6.3, 0.5 },  { 5, -0.4, 1.4, -2 },  { 8.4, 5, 7.4, 2.3 } }, "DNA", true, "Packer, M. J. ; Dauncey, M. P. ; Hunter, C. A. (2000): Sequence-dependent DNA Structure: Dinucleotide Conformational Maps.  J. Mol. Biol. (2000) 295 (1), 71-83", "10623509", HowCreated.CALCULATED, Type.CONFORMATIONAL, "Dimension", "Optimised values of roll at zero slide and shift." ),

	TWIST_PACKER(new double[][]{ { 37.6, 35.8, 35.7, 39.7 },  { 32.2, 35.5, 33.9, 35.7 },  { 38.4, 37.4, 35.5, 35.8 },  { 34.6, 38.4, 32.2, 37.6 } }, "DNA", true, "Packer, M. J. ; Dauncey, M. P. ; Hunter, C. A. (2000): Sequence-dependent DNA Structure: Dinucleotide Conformational Maps.  J. Mol. Biol. (2000) 295 (1), 71-83", "10623509", HowCreated.CALCULATED, Type.CONFORMATIONAL, "Dimension", "Optimised values of twist at zero slide and shift." ),

	FLEXIBILITY_SLIDE(new double[][]{ { 13.72, 9.57, 7.58, 11.69 },  { 1.35, 7.36, 4.02, 7.58 },  { 10.28, 4.34, 7.36, 9.57 },  { 7.13, 10.28, 1.35, 13.72 } }, "DNA", true, "Packer, M. J. ; Dauncey, M. P. ; Hunter, C. A. (2000): Sequence-dependent DNA Structure: Dinucleotide Conformational Maps.  J. Mol. Biol. (2000) 295 (1), 71-83  ", "10623509", HowCreated.CALCULATED, Type.PHYSICOCHEMICAL, "kJ mol^-1 A^-2", "Flexibility F_slide determined from the curvature of the slide/shift stacking potential at the minimum energy. This amounts to a force constant for sliding the step." ),

	FLEXIBILITY_SHIFT(new double[][]{ { 5.35, 9.73, 8.98, 1.13 },  { 4.61, 5.51, 12.13, 8.98 },  { 5.44, 1.98, 5.51, 9.73 },  { 4.28, 5.44, 4.61, 5.35 } }, "DNA", true, "Packer, M. J. ; Dauncey, M. P. ; Hunter, C. A. (2000): Sequence-dependent DNA Structure: Dinucleotide Conformational Maps.  J. Mol. Biol. (2000) 295 (1), 71-83", "10623509", HowCreated.CALCULATED, Type.PHYSICOCHEMICAL, "kJ mol^-1 A^-2", "Flexibility F_shift determined from the curvature of the slide/shift stacking potential at the minimum energy. This amounts to a force constant for shifting the step." ),

	ENTHALPY_SANTALUCIA(new double[][]{ { -7.6, -8.4, -7.8, -7.2 },  { -8.5, -8, -10.6, -7.8 },  { -8.2, -9.8, -8, -8.4 },  { -7.2, -8.2, -8.5, -7.6 } }, "DNA", true, "SantaLucia, J.  ;  Hicks, D. (2004): The Thermodynamics of DNA Structural Motifs  Annu. Rev. Biophys. Biomol. Struct. (2004) 33, 415?40", "15139820", HowCreated.CALCULATED, Type.PHYSICOCHEMICAL, "kcal/mol", "delta H" ),

	ENTROPY_SANTALUCIA(new double[][]{ { -21.3, -22.4, -21, -20.4 },  { -22.7, -19.9, -27.2, -21 },  { -22.2, -24.4, -19.9, -22.4 },  { -21.3, -22.2, -22.7, -21.3 } }, "DNA", true, "SantaLucia, J.  ;  Hicks, D. (2004): The Thermodynamics of DNA Structural Motifs  Annu. Rev. Biophys. Biomol. Struct. (2004) 33, 415?40", "15139820", HowCreated.CALCULATED, Type.PHYSICOCHEMICAL, "e.u.", "delta S" ),

	FREE_ENERGY_SANTALUCIA_2004(new double[][]{ { -1, -1.44, -1.28, -0.88 },  { -1.45, -1.84, -2.17, -1.28 },  { -1.3, -2.24, -1.84, -1.44 },  { -0.58, -1.3, -1.45, -1 } }, "DNA", true, "SantaLucia, J.  ;  Hicks, D. (2004): The Thermodynamics of DNA Structural Motifs  Annu. Rev. Biophys. Biomol. Struct. (2004) 33, 415?40", "15139820", HowCreated.CALCULATED, Type.PHYSICOCHEMICAL, "kcal/mol", "delta G" );
	
	public static final Smoothing NO_SMOOTHING = new NoSmoothing();
	
	public static final AlphabetContainer continuousAlphabet = new AlphabetContainer( new ContinuousAlphabet() );

	/**
	 * This enum defines the types of dinucleotide properties.
	 * @author Jan Grau
	 *
	 */
	public enum Type{
		/**
		 * A conformational property of nucleotide sequences
		 */
		CONFORMATIONAL,
		/**
		 * A physicochemical property of nucleotide sequences
		 */
		PHYSICOCHEMICAL,
		/**
		 * A letter-based property of nucleotide sequences
		 */
		LETTER_BASED;
	}
	
	/**
	 * This enum defines the origins of nucleotide properties
	 * @author Jan Grau
	 *
	 */
	public enum HowCreated{
		/**
		 * The property has been determined experimentally
		 */
		EXPERIMENTAL,
		/**
		 * The property has been calculated
		 */
		CALCULATED;
	}

	
	private enum Annotation{
		NO_ANNOTATION,
		ORIGINAL_AS_ANNOTATION,
		ADD_PROPERTY_AS_ANNOTATION,
		SET_PROPERTY_AS_ANNOTATION;
	}
	
	private double[][] dinucleotideParameters;
	private String nucleicAcid;
	private boolean doubleStrand;
	private String reference;
	private String pubMedID;
	private HowCreated howCreated;
	private Type type;
	private String dimension;
	private String comments;
	
	
	private DinucleotideProperty(double[][] dinucleotideParameters, String nucleicAcid, boolean doubleStrand, String reference, String pubMedID, HowCreated howCreated, Type type, String dimension, String comments){
		this.dinucleotideParameters = dinucleotideParameters;
		this.nucleicAcid = nucleicAcid;
		this.doubleStrand = doubleStrand;
		this.reference = reference;
		this.pubMedID = pubMedID;
		this.howCreated = howCreated;
		this.type = type;
		this.dimension = dimension;
		this.comments = comments;
	}

	/**
	 * Computes this dinucleotide property for all overlapping twomers in <code>original</code>, smoothes the result using <code>smoothing</code>,
	 * and returns the smoothed property as a <code>double</code> array. The length of this <code>double</code> array depends on the smoothing applied.
	 * @param original the original nucleotide sequence
	 * @param smoothing the smoothing applied to the property values
	 * @return the smoothed property for all overlapping twomers
	 * @throws WrongSequenceTypeException if <code>original</code> is not a DNA sequence
	 */
	public double[] getProperty(Sequence original, Smoothing smoothing) throws WrongSequenceTypeException{
		double[] prop = new double[original.getLength()-1];
		return getProperty( original, prop, smoothing );
	}
		
	private double[] getProperty(Sequence original, double[] prop, Smoothing smoothing) throws WrongSequenceTypeException{
		if(!original.getAlphabetContainer().isSimple() || !original.getAlphabetContainer().getAlphabetAt( 0 ).checkConsistency( DNAAlphabet.SINGLETON )){
			throw new WrongSequenceTypeException();
		}
		int curr;
		int last = original.discreteVal( 0 );
		for(int i=0;i<prop.length;i++){
			curr = original.discreteVal( i+1 );
			prop[i] = dinucleotideParameters[last][curr];
			last = curr;
		}
		
		prop = smoothing.smooth( prop );
		return prop;
	}
	
	public double getProperty( Sequence seq, int start ) {
		if(seq instanceof SubSequence){
			start += ((SubSequence) seq).getStart();
			seq = ((SubSequence) seq).getOriginal();
		}
		if(start == 0){
			double mean = 0;
			for(int i=0;i<dinucleotideParameters.length;i++){
				mean += dinucleotideParameters[i][seq.discreteVal( start )];
			}
			mean /= dinucleotideParameters.length;
			return mean;
		}else{
			return dinucleotideParameters[seq.discreteVal( start-1 )][seq.discreteVal( start )];
		}
	}
	
	/**
	 * Computes this dinucleotide property for all overlapping twomers in <code>original</code>
	 * and returns the result as a <code>double</code> array of length <code>original.getLength()-1</code>
	 * @param original the original nucleotide sequence
	 * @return the property for all overlapping twomers
	 * @throws WrongSequenceTypeException if <code>original</code> is not a DNA sequence
	 */
	public double[] getProperty(Sequence original) throws WrongSequenceTypeException{
		return getProperty( original, NO_SMOOTHING );
	}
	
	/**
	 * Computes this dinucleotide property for all overlapping dimers in <code>original</code>
	 * and returns the result as a {@link Sequence} of length <code>original.getLength()-1</code>
	 * @param original the original nucleotide sequence
	 * @return the property for all overlapping dimers
	 * @throws WrongAlphabetException if the new sequence could not be created
	 * @throws WrongSequenceTypeException if <code>original</code> is not a DNA sequence
	 */
	public ArbitrarySequence getPropertyAsSequence(Sequence original) throws WrongAlphabetException, WrongSequenceTypeException{
		return getPropertyAsSequence( original, NO_SMOOTHING );
	}
	
	/**
	 * Computes this dinucleotide property for all overlapping dimers in <code>original</code>, smoothes the result using <code>smoothing</code>,
	 * and returns the smoothed property as a {@link Sequence}. The length of this {@link Sequence} depends on the smoothing applied.
	 * @param original the original nucleotide sequence
	 * @param smoothing the smoothing applied to the property values
	 * @return the smoothed property for all overlapping dimers
	 * @throws WrongSequenceTypeException if <code>original</code> is not a DNA sequence
	 */
	public ArbitrarySequence getPropertyAsSequence(Sequence original, Smoothing smoothing) throws WrongSequenceTypeException{
		try{
			return new ArbitrarySequence( continuousAlphabet, getProperty( original, smoothing ) );
		}catch(WrongAlphabetException doesnothappen){
			return null;
		}
		
	}
	
	/**
	 * Returns the dinucleotide parameters of this {@link DinucleotideProperty} as a two-dimensional <code>double</code> array, where the rows correspond to the first nucleotide
	 * and the columns correspond to the second nucleotide in the dinucleotide in order A, C, G, and T.
	 * @return the dinucleotide parameters
	 */
	public double[][] getDinucleotideParameters() {
		return dinucleotideParameters.clone();
	}


	/**
	 * Returns the kind of nucleic acid, e.g. RNA, DNA, this property has been determined for.
	 * @return the kind of nucleic acid
	 */
	public String getNucleicAcid() {
		return nucleicAcid;
	}


	/**
	 * Returns <code>true</code> if this property has been determined for a double-stranded nucleic acid.
	 * @return <code>true</code> if this property has been determined for a double-stranded nucleic acid
	 */
	public boolean isDoubleStrand() {
		return doubleStrand;
	}


	/**
	 * Returns the reference of the publication where the parameters of this property has been published.
	 * @return the reference
	 */
	public String getReference() {
		return reference;
	}


	/**
	 * Returns the PubMed ID of the publication where the parameters of this property has been published.
	 * @return the PubMed ID
	 */
	public String getPubMedID() {
		return pubMedID;
	}

	
	/**
	 * Returns how this property has been determined.
	 * @return one of {@link HowCreated}
	 */
	public HowCreated getHowCreated() {
		return howCreated;
	}


	/**
	 * Returns the type of this property.
	 * @return one of {@link Type}
	 */
	public Type getType() {
		return type;
	}


	/**
	 * Returns the dimension of this property, e.g. kcal/mol or angstroem.
	 * @return the dimension
	 */
	public String getDimension() {
		return dimension;
	}


	/**
	 * Returns additional comments on this property.
	 * @return the comments
	 */
	public String getComments() {
		return comments;
	}
	
	
	/**
	 * Filters all {@link DinucleotideProperty}s by some of their annotations.
	 * @param nucleicAcid the kind of the nucleic acid or <code>null</code> if no filter shall be applied
	 * @param doubleStrand the strandedness or <code>null</code> if no filter shall be applied
	 * @param howCreated the kind of creation or <code>null</code> if no filter shall be applied
	 * @param type the type of the property or <code>null</code> if no filter shall be applied
	 * @param dimension the dimension of the property or <code>null</code> if no filter shall be applied
	 * @return a filtered subset of {@link DinucleotideProperty}s
	 */
	public static DinucleotideProperty[] filterProperties(String nucleicAcid, Boolean doubleStrand, HowCreated howCreated, Type type, String dimension){
		DinucleotideProperty[] vals = DinucleotideProperty.values();
		
		int num = vals.length;
		
		for(int i=0;i<vals.length;i++){
			if(nucleicAcid != null && !vals[i].getNucleicAcid().equals( nucleicAcid )){
				vals[i] = null;
				num--;
			}
			if(doubleStrand != null && vals[i] != null && vals[i].isDoubleStrand() != doubleStrand){
				vals[i] = null;
				num--;
			}
			if(howCreated != null && vals[i] != null && vals[i].getHowCreated() != howCreated){
				vals[i] = null;
				num--;
			}
			if(type != null && vals[i] != null && vals[i].getType() != type){
				vals[i] = null;
				num--;
			}
			if(dimension != null && vals[i] != null && !vals[i].getDimension().equals( dimension )){
				vals[i] = null;
				num--;
			}
		}
		
		DinucleotideProperty[] sel = new DinucleotideProperty[num];
		
		for(DinucleotideProperty curr : vals){
			if(curr != null){
				sel[sel.length-num] = curr;
				num--;
			}
		}
		
		return sel;
	}
	
	/**
	 * Creates a new {@link DataSet} by converting each {@link Sequence} in <code>original</code> to the {@link DinucleotideProperty} <code>property</code>.
	 * @param original {@link DataSet} containing the original {@link Sequence}s
	 * @param property the property
	 * @return the converted {@link Sequence}s as a {@link DataSet}
	 * @throws WrongSequenceTypeException if <code>original</code> contains non-DNA sequences
	 */
	public static DataSet getDataSetForProperty(DataSet original, DinucleotideProperty property) throws WrongSequenceTypeException{
		return getDataSetForProperty( original, NO_SMOOTHING, false, property );
	}
	
	/**
	 * Creates a new {@link DataSet} by converting each {@link Sequence} in <code>original</code> to the {@link DinucleotideProperty} <code>property</code> using the {@link Smoothing} smoothing.
	 * @param original {@link DataSet} containing the original {@link Sequence}s
	 * @param smoothing the smoothing
	 * @param originalAsAnnotation if <code>true</code>, the original {@link Sequence} is added as a {@link ReferenceSequenceAnnotation} to each converted {@link Sequence}
	 * @param property the property
	 * @return the converted {@link Sequence}s as a {@link DataSet}
	 * @throws WrongSequenceTypeException if <code>original</code> contains non-DNA sequences
	 */
	public static DataSet getDataSetForProperty(DataSet original, Smoothing smoothing, boolean originalAsAnnotation, DinucleotideProperty property) throws WrongSequenceTypeException{
		return getDataSetForProperty( original, smoothing, originalAsAnnotation ? Annotation.ORIGINAL_AS_ANNOTATION : Annotation.NO_ANNOTATION, property );
	}
	
	/**
	 * Creates a new {@link DataSet} by converting each {@link Sequence} in <code>original</code> to the {@link DinucleotideProperty}s <code>properties</code> and setting these as {@link ReferenceSequenceAnnotation} of each original sequence.
	 * @param original {@link DataSet} containing the original {@link Sequence}s
	 * @param properties the properties
	 * @return the annotated {@link Sequence}s as a {@link DataSet}
	 * @throws WrongSequenceTypeException if <code>original</code> contains non-DNA sequences
	 */
	public static DataSet getDataSetForProperty(DataSet original, DinucleotideProperty... properties) throws WrongSequenceTypeException{
		return getDataSetForProperty( original, NO_SMOOTHING, Annotation.SET_PROPERTY_AS_ANNOTATION, properties );
	}
	
	/**
	 * Creates a new {@link DataSet} by converting each {@link Sequence} in <code>original</code> to the {@link DinucleotideProperty}s <code>properties</code> and adding or setting these as {@link ReferenceSequenceAnnotation} of each original sequence.
	 * @param original {@link DataSet} containing the original {@link Sequence}s
	 * @param smoothing the smoothing
	 * @param addToAnnotation if <code>true</code> the converted {@link Sequence}s are added to the current annotation, otherwise the current annotation is replaced
	 * @param properties the properties
	 * @return the {@link DataSet} of property sequences
	 * @throws WrongSequenceTypeException if <code>original</code> contains non-DNA sequences
	 */
	public static DataSet getDataSetForProperty(DataSet original, Smoothing smoothing, boolean addToAnnotation, DinucleotideProperty... properties) throws WrongSequenceTypeException{
		return getDataSetForProperty( original, smoothing, addToAnnotation ? Annotation.ADD_PROPERTY_AS_ANNOTATION : Annotation.SET_PROPERTY_AS_ANNOTATION, properties );
	}
	
	private static DataSet getDataSetForProperty(DataSet original, Smoothing smoothing, Annotation annotation, DinucleotideProperty... properties) throws WrongSequenceTypeException{
		Sequence[] seqs = new Sequence[original.getNumberOfElements()];
		
		if(properties.length > 1 && annotation == Annotation.NO_ANNOTATION || annotation == Annotation.ORIGINAL_AS_ANNOTATION){
			throw new UnsupportedOperationException();
		}
		
		StringBuffer ann = new StringBuffer();
		String[] anns = new String[properties.length];
		for(int i=0;i<anns.length;i++){
			ann.append(properties[i].name());
			if(i < anns.length-1){
				ann.append( ", " );
			}
			anns[i] = properties[i].name()+" with smoothing "+smoothing.toString();
		}
		ann.append( " with smoothing "+smoothing.toString() );
		
		for(int i=0;i<seqs.length;i++){
			Sequence seq = original.getElementAt( i );
			if(annotation == Annotation.NO_ANNOTATION || annotation == Annotation.ORIGINAL_AS_ANNOTATION){
				seqs[i] = properties[0].getPropertyAsSequence( seq, smoothing );
				if(annotation == Annotation.ORIGINAL_AS_ANNOTATION){
					seqs[i] = seqs[i].annotate( false, new ReferenceSequenceAnnotation( "original", seq ) );
				}
			}else{
				SequenceAnnotation[] props = new SequenceAnnotation[properties.length];
				for(int j=0;j<props.length;j++){
					props[j] = new ReferenceSequenceAnnotation(anns[j], properties[j].getPropertyAsSequence( seq, smoothing ));
				}
				seqs[i] = seq.annotate( annotation == Annotation.ADD_PROPERTY_AS_ANNOTATION, props );
				
			}
		}
		if(annotation == Annotation.NO_ANNOTATION || annotation == Annotation.ORIGINAL_AS_ANNOTATION){
			try{
				return new DataSet( original.getAnnotation()+" converted to "+ann.toString(), seqs );
			}catch(EmptyDataSetException doesnothappen){
				return null;
			}catch(WrongAlphabetException doesnothappen2){
				return null;
			}
		}else{
			try{
				return new DataSet( original.getAnnotation()+" with annotation "+ann.toString(), seqs );
			}catch(EmptyDataSetException doesnothappen){
				return null;
			}catch(WrongAlphabetException doesnothappen2){
				return null;
			}
			
		}
	}
	
	public static ImageResult getPropertyImage( Sequence original, DinucleotideProperty prop, Smoothing smoothing, REnvironment re, int xLeft, String pltOptions, int width, int height ) throws Exception {
		double[] values = prop.getProperty( original, smoothing );
		re.createVector( "values", values );
		if( pltOptions == null ) {
			pltOptions = "\"l\",xlab=\"position\",ylab=\""+prop.getDimension()+"\",las=1,main=\"" + prop.name() + "\"";
		}
		BufferedImage img = re.plot( "plot(1:length(values)-1+"+xLeft+",values," + pltOptions + ");", width, height );
		return new ImageResult("profile for " + prop.name(), "the profile for the dinucleotide property " + prop.name() + " and the smoothing " + smoothing.toString(), img );
	}
	
	public static ImageResult getPropertyImage( DataSet original, DinucleotideProperty prop, Smoothing smoothing, REnvironment re, int xLeft, String pltOptions, int width, int height ) throws Exception {
		int l = original.getElementLength();
		if( l == 0 ) {
			throw new WrongLengthException( "All Sequences of the DataSet have to have the same length." );
		} else {
			l--; 
		}
		if( pltOptions == null ) {
			pltOptions = "ylim=c(min(values),max(values)),xlab=\"position\",ylab=\""+prop.getDimension()+"\",las=1,main=\"" + prop.name() + "\"";
		}
		double[] help = new double[l];
		double[][] values = new double[original.getNumberOfElements()][];
		int i = 0;
		for( Sequence seq : original ) {
			values[i] = prop.getProperty( seq, help, smoothing );
			if( values[i] == help ) {
				values[i] = values[i].clone();
			}
			i++;
		}
		re.createMatrix( "values", values );
		re.voidEval( "d=dim(values);" );
		re.voidEval( "x=1:(d[2])-1+"+xLeft );
		
		String pltcmd = "plot(x,values[1,],col=0," + pltOptions + ");\n"
				+ "n=d[1];for(i in 1:n){lines(x,values[i,],col=gray(0.8),lty=2);}\n"
				+ "lines(x,apply(values,2,mean),col=2);"
		;
		BufferedImage img = re.plot( pltcmd, width, height );
		return new ImageResult("profile for " + prop.name(), "the profile for the dinucleotide property " + prop.name() + " and the smoothing " + smoothing.toString(), img );
	}

/*
	public static void filter() {
		DinucleotideProperty[] all = values();
		int anz = 0;
		for( DinucleotideProperty current : all ) {
			int[] b1 = current.check( 1 );
			int[] b2 = current.check( -1 );
			if( ! (b1==null || b2==null) ) {
				System.out.println( current + "\t" + current.doubleStrand + "\t"
						+ "\t" + Arrays.toString( current.dinucleotideParameters[0] )
						+ "\t" + Arrays.toString( current.dinucleotideParameters[1] )
						+ "\t" + Arrays.toString( current.dinucleotideParameters[2] )
						+ "\t" + Arrays.toString( current.dinucleotideParameters[3] ) );
				if( b1 != null ) {
					System.out.println( "1\t" + Arrays.toString( b1 ) + "\t"
							+ current.dinucleotideParameters[b1[0]][b1[1]] + "\t" + current.dinucleotideParameters[3-b1[1]][3-b1[0]] );
				}
				if( b2 != null ) {
					System.out.println( "-1\t" + Arrays.toString( b2 ) + "\t"
							+ current.dinucleotideParameters[b2[0]][b2[1]] + "\t" + current.dinucleotideParameters[3-b2[1]][3-b2[0]] );
				}
				anz++;
			}
		}
		System.out.println( anz );
	}
	
	private int[] check( double factor ) {
		int all = 3, outer = -1, inner = all+1;
		//System.out.println(name() + "\t" + factor );
		while( outer < all && inner > all ) {
			outer++;
			inner = 0;
			while( inner <= all && dinucleotideParameters[outer][inner] == factor*dinucleotideParameters[all-inner][all-outer] ) {
				//System.out.println( "(" + outer + "," + inner+ ")\t" + "(" + (all-inner) + "," + (all-outer)+ ")\t" + dinucleotideParameters[outer][inner] + " == " + (factor*dinucleotideParameters[all-inner][all-outer]) + "\t" + (dinucleotideParameters[outer][inner] == factor*dinucleotideParameters[all-inner][all-outer]) );
				inner++;
			}
			//System.out.println( outer + "\t" + inner );
		}
		//System.out.println();
		return outer == all ? null : new int[]{outer, inner};
	}
/**/	
	
	/**
	 * Abstract class for methods that smooth a series of real values.
	 * @author Jan Grau
	 *
	 */
	public static abstract class Smoothing{
		
		/**
		 * Returns the smoothed version of <code>original</code>. The length of the returned array and <code>original</code> may differ.
		 * @param original the original series of values
		 * @return the smoothed version
		 */
		public abstract double[] smooth(double[] original);
		
		/**
		 * Returns a {@link String} representation of this smoothing method.
		 */
		public abstract String toString();
		
	}
	
	/**
	 * Implementation of {@link Smoothing} that conducts no smoothing. The static instance {@link DinucleotideProperty#NO_SMOOTHING} should be used instead of constructing new objects.
	 * @author Jan Grau
	 *
	 */
	public static class NoSmoothing extends Smoothing{
		
		private NoSmoothing(){}
		
		/* (non-Javadoc)
		 * @see de.jstacs.data.DinucleotideProperty.Smoothing#smooth(double[])
		 */
		public double[] smooth(double[] original){
			return original;
		}
		
		/* (non-Javadoc)
		 * @see de.jstacs.data.DinucleotideProperty.Smoothing#toString()
		 */
		public String toString(){
			return "none";
		}
		
	}
	
	
	/**
	 * Smoothing by mean using a pre-defined window width.
	 * 
	 * @author Jan Grau
	 *
	 */
	public static class MeanSmoothing extends Smoothing{
		
		private int width;
		
		/**
		 * Creates a new {@link MeanSmoothing} that averages over windows of width <code>width</code>.
		 * @param width the width of the window.
		 */
		public MeanSmoothing(int width){
			this.width = width;
		}
		
		/* (non-Javadoc)
		 * @see de.jstacs.data.DinucleotideProperty.Smoothing#smooth(double[])
		 */
		public double[] smooth(double[] original){
			if(original.length < width){
				return new double[0];
			}
			double[] smoothed = new double[original.length-width+1];
			
			double curr = 0;
			int i = 0, w = width-1;
			for(i=0;i<w;i++){
				curr += original[i];
			}
			
			for(i=0;i<smoothed.length;i++){
				curr += original[i+w];
				smoothed[i] = curr/width;
				curr -= original[i];
			}
			return smoothed;
		}

		/* (non-Javadoc)
		 * @see de.jstacs.data.DinucleotideProperty.Smoothing#toString()
		 */
		@Override
		public String toString() {
			return "mean smoothing with width "+width;
		}
		
	}
	

	/**
	 * Smoothing by median using a pre-defined window width.
	 * 
	 * @author Jan Grau
	 *
	 */
	public static class MedianSmoothing extends Smoothing{
		
		private int width;
		private double[] sorted;
		
		/**
		 * Creates a new {@link MedianSmoothing} that computes the median over windows of width <code>width</code>.
		 * @param width the width of the window.
		 */
		public MedianSmoothing(int width){
			this.width = width;
			this.sorted = new double[width];
		}
		
		/* (non-Javadoc)
		 * @see de.jstacs.data.DinucleotideProperty.Smoothing#smooth(double[])
		 */
		public double[] smooth(double[] original){
			if(original.length < width){
				return new double[0];
			}
			double[] smoothed = new double[original.length-width+1];
			
			System.arraycopy( original, 0, sorted, 0, width );
			Arrays.sort( sorted );
			
			int idx = -1, w = width-1;
			double old = original[w];
			
			if(width % 2 == 0){
				idx = width/2;
				for(int i=0;i<smoothed.length;i++){
					bubble(sorted,old,original[i+w]);
					old = original[i];
					smoothed[i] = (sorted[idx-1] + sorted[idx])/2.0;
				}
			}else{
				idx = (width-1)/2;
				for(int i=0;i<smoothed.length;i++){
					bubble(sorted,old,original[i+w]);
					old = original[i];
					smoothed[i] = sorted[idx];					
				}
			}
			return smoothed;
		}
		
		private static void bubble(double[] vals, double oldVal, double newVal){
			if(oldVal == newVal){
				return;
			}
			int idx =  Arrays.binarySearch( vals, oldVal );
			if(vals[idx] > newVal){
				while(idx > 0 && vals[idx-1] > newVal){
					vals[idx] = vals[idx-1];
					idx--;
				}
				vals[idx] = newVal;
			}else{
				while(idx < vals.length-1 && vals[idx+1] < newVal){
					vals[idx] = vals[idx+1];
					idx++;
				}
				vals[idx] = newVal;
			}
		}
		
		/* (non-Javadoc)
		 * @see de.jstacs.data.DinucleotideProperty.Smoothing#toString()
		 */
		@Override
		public String toString() {
			return "median smoothing with width "+width;
		}
		
	}
	
}
