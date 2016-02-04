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

package projects.slim;

import java.io.File;
import java.util.LinkedList;

import javax.imageio.ImageIO;

import projects.dimont.Interpolation;
import projects.dimont.ThresholdedStrandChIPper;
import de.jstacs.classifiers.differentiableSequenceScoreBased.gendismix.GenDisMixClassifier;
import de.jstacs.data.DataSet;
import de.jstacs.data.alphabets.DNAAlphabet;
import de.jstacs.data.alphabets.DNAAlphabetContainer;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.data.sequences.SparseSequence;
import de.jstacs.data.sequences.annotation.SequenceAnnotation;
import de.jstacs.data.sequences.annotation.SequenceAnnotationParser;
import de.jstacs.data.sequences.annotation.SplitSequenceAnnotationParser;
import de.jstacs.io.FileManager;
import de.jstacs.motifDiscovery.SignificantMotifOccurrencesFinder;
import de.jstacs.parameters.ParameterSetTagger;
import de.jstacs.results.ImageResult;
import de.jstacs.results.ListResult;
import de.jstacs.results.Result;
import de.jstacs.utils.DoubleList;
import de.jstacs.utils.PFMComparator;
import de.jstacs.utils.Pair;
import de.jstacs.utils.SeqLogoPlotter;
import de.jstacs.utils.ToolBox;


public class SlimDimontPredictor {

	
	
	/**
	 * @param args
	 */
	public static void main( String[] args ) throws Exception {
		
		ParameterSetTagger cParams = new ParameterSetTagger( SlimDimontPredictorParameterSet.PREFIX, new SlimDimontPredictorParameterSet() );
		cParams.fillParameters( "=", args );
		System.out.println( "parameters:" );
		System.out.println( cParams );
		System.out.println("_________________________________");
		if( !cParams.hasDefaultOrIsSet() ) {
			System.out.println( "Some of the required parameters are not specified." );
			System.exit( 1 );
		}

		
		
		
		String home = cParams.getValueFromTag( SlimDimontPredictorParameterSet.HOME, String.class );
		String fgData = home +File.separator+ cParams.getValueFromTag( SlimDimontPredictorParameterSet.DATA, String.class );
		String infix = cParams.getValueFromTag( SlimDimontPredictorParameterSet.INFIX, String.class );
		String value = cParams.getValueFromTag( SlimDimontParameterSet.VALUE_TAG, String.class );
		String weightingFactor = cParams.getValueFromTag( SlimDimontParameterSet.WEIGHTING_FACTOR, String.class );
		
		String chipper = cParams.getValueFromTag( SlimDimontPredictorParameterSet.CHIPPER, String.class );
		double pval = cParams.getValueFromTag( SlimDimontPredictorParameterSet.PVAL, Double.class );
		
		GenDisMixClassifier cl = new GenDisMixClassifier( FileManager.readFile( chipper ) );
		
		ThresholdedStrandChIPper model = (ThresholdedStrandChIPper)cl.getDifferentiableSequenceScore( 0 );
		
		SequenceAnnotationParser parser = new SplitSequenceAnnotationParser(":", ";");
		
		Result[][] res = run(SparseSequence.getDataSet(DNAAlphabetContainer.SINGLETON, fgData, parser),value,weightingFactor,model,pval);
		
		for(int i=0;i<res.length;i++){
			System.out.println("+++++++++++++++++++++++++++++++++++++++++++++\nMotif model:");
			System.out.println( model);
			System.out.println("+++++++++++++++++++++++++++++++++++++++++++++\n");

			ListResult lr = (ListResult) res[i][0];
			System.out.println("Predicted "+lr.getRawResult().length+" binding sites");
			FileManager.writeFile( new File(home+File.separator+infix+"-predictions.txt"), lr.toString() );
			
			if(res[i].length > 1){
				ImageResult ir = (ImageResult) res[i][1];
				ImageIO.write( ir.getValue() , "png", new File(home+File.separator+infix+"-logo.png") );
				ir = (ImageResult) res[i][2];
				ImageIO.write( ir.getValue() , "png", new File(home+File.separator+infix+"-logo-rc.png") );
				ir = (ImageResult) res[i][3];
				ImageIO.write( ir.getValue() , "png", new File(home+File.separator+infix+"-dependencylogo.png") );
			}
			
		}
		
	}

	public static Result[][] run( DataSet data, String value, String weightingFactor, ThresholdedStrandChIPper model, double pval ) throws Exception {
		
		double[][] weights = new double[2][data.getNumberOfElements()];
		
		double[] raw = weights[0].clone();
		
		//read annotation
		for( int j = 0; j < weights[0].length; j++ ) {
			Sequence seq = data.getElementAt(j);
			SequenceAnnotation[] seqAn = seq.getAnnotation();
			for( int i = 0; i < seqAn.length; i++ ) {
				if( seqAn[i].getType().equals(value) ) {
					raw[j] = Double.parseDouble( seqAn[i].getIdentifier() );
				}
			}
		}
		
		//create weights
		double wf;
		if( weightingFactor.endsWith("sd") ) {
			double h = Double.parseDouble( weightingFactor.substring(0,weightingFactor.length()-2) );
			double meanRaw = ToolBox.sum(raw) / raw.length;
			double sdRaw = 0;
			for( int i = 0; i < raw.length; i++ ) {
				sdRaw += (raw[i]-meanRaw) * (raw[i]-meanRaw);
			}
			sdRaw = Math.sqrt( sdRaw/raw.length );
			h = meanRaw + h*sdRaw;
			double anz = 0;
			for( int i = 0; i < raw.length; i++ ) {
				if( raw[i] >= h ) {
					anz++;
				}
			}
			anz=Math.max(50,anz);
			wf = anz/raw.length;
		} else {
			wf = Double.parseDouble( weightingFactor );
		}

		weights[0] = Interpolation.getWeight( data, raw, wf, Interpolation.RANK_LOG );
		weights[1] = Interpolation.getBgWeight( weights[0] );
		
	
		
		SignificantMotifOccurrencesFinder smof = new SignificantMotifOccurrencesFinder( model, data, weights[1], pval );
		
		Pair<double[][][],int[][]> pair = smof.getPWMAndPositions( 0, data, weights[0], 0, 0 );
		
		LinkedList<Result> result = new LinkedList<Result>();
		
		
		LinkedList<Sequence> bs = new LinkedList<Sequence>();
		DoubleList bsWeights = new DoubleList();
		
		result.add(SlimDimont.getListResult(data, weights[0], pair, model.getMotifLength( 0 ), 0, model.getMotifModel(), bs, bsWeights, value ));
		
		double[][] pwm = pair.getFirstElement()[0];
		
		if(!Double.isNaN( pwm[0][0] )){
			try{
				int height = SeqLogoPlotter.getHeight( 750, pwm );
				result.add(new ImageResult( "Motif", "Sequence logo of the motif", SeqLogoPlotter.plotLogoToBufferedImage( height, pwm ) ));
				result.add(new ImageResult( "Motif (rc)", "Sequence logo of the reverse complement of the motif", SeqLogoPlotter.plotLogoToBufferedImage( height, PFMComparator.getReverseComplement( DNAAlphabet.SINGLETON, pwm ) ) ));
				result.add( new ImageResult( "Dependency logo", "Dependency logo of the motif", SeqLogoPlotter.plotDefaultDependencyLogoToBufferedImage( new DataSet("",bs), bsWeights.toArray(), 600 ) ) );
			}catch(Exception e){
				
			}catch(InternalError er){
				
			}
		}
		
		return new Result[][]{result.toArray( new Result[0] )};
		
	}

}
