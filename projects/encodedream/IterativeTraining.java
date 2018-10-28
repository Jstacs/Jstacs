package projects.encodedream;
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
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;

import de.jstacs.algorithms.optimization.Optimizer;
import de.jstacs.classifiers.differentiableSequenceScoreBased.OptimizableFunction.KindOfParameter;
import de.jstacs.classifiers.differentiableSequenceScoreBased.gendismix.GenDisMixClassifier;
import de.jstacs.classifiers.differentiableSequenceScoreBased.gendismix.GenDisMixClassifierParameterSet;
import de.jstacs.classifiers.differentiableSequenceScoreBased.gendismix.LearningPrinciple;
import de.jstacs.classifiers.differentiableSequenceScoreBased.logPrior.DoesNothingLogPrior;
import de.jstacs.data.DataSet;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.sequenceScores.statisticalModels.differentiable.continuous.GaussianNetwork;
import de.jstacs.utils.DoubleList;
import de.jstacs.utils.Pair;
import de.jstacs.utils.ToolBox;

public class IterativeTraining {
	
	private FeatureReader reader;
	private int threads;
	private HashMap<String,Integer> sizes;
	
	public IterativeTraining(FeatureReader reader, int threads, HashMap<String,Integer> sizes){
		this.reader = reader;
		this.threads = threads;
		this.sizes = sizes;
	}
	
	
	public GenDisMixClassifier train(DataSet[] data, double[][] weights) throws Exception{
		
		GaussianNetwork gn = new GaussianNetwork(new int[data[0].getElementLength()][0]);
		
		GenDisMixClassifierParameterSet params = new GenDisMixClassifierParameterSet(data[0].getAlphabetContainer(), gn.getLength(), Optimizer.QUASI_NEWTON_BFGS, 1E-6, 1E-6, 1E-4, false, KindOfParameter.PLUGIN, true, threads);
		
		GenDisMixClassifier cl = new GenDisMixClassifier(params, DoesNothingLogPrior.defaultInstance, LearningPrinciple.MCL, gn,gn);
		
		cl.train(data, weights);
		
		return cl;
		
	}
	
	
	
	public GenDisMixClassifier[] iterativeTraining(int iterations, HashSet<String> trainChroms, LinkedList<String> itChroms, double perc, int binsBefore, int binsAfter) throws Exception{
		if(trainChroms!= null && !trainChroms.containsAll(itChroms)){
			throw new RuntimeException();
		}
		
		Pair<DataSet[], double[][]> pair = reader.getInitialData(trainChroms);
		
		
		LinkedList<GenDisMixClassifier> clList = new LinkedList<>();
		GenDisMixClassifier cl = train(pair.getFirstElement(),pair.getSecondElement());
		clList.add(cl);
		
		DataSet fg = pair.getFirstElement()[0];
		DataSet bg = pair.getFirstElement()[1];
		double[][] weights = pair.getSecondElement();
		
		for(int i=1;i<iterations;i++){
			
//System.out.println("Iteration "+i);
			Predictor pred = new Predictor(clList.toArray(new GenDisMixClassifier[0]), reader,binsBefore,binsAfter);
			//double[] scPos = pred.predict(fg);
//System.out.println("Predicted");
			//double th = ToolBox.percentile(scPos, perc);
//System.out.println("threshold: "+th);
			
			DoubleList scPos = new DoubleList();
			DoubleList scNeg = new DoubleList();
			
			double[][] preds = new double[itChroms.size()][];
			
			for(int l=0;l<itChroms.size();l++){
				String chr = itChroms.get(l);
				preds[l] = pred.predict(chr, sizes.get(chr));
			}
			
			
			
			reader.reset();
			for(int l=0;l<itChroms.size();l++){
				String chr = itChroms.get(l);
				reader.findChr(chr);
				int j=0;
				int size = sizes.get(chr);
				do{
					if(j<size){
						
						
						char lab = reader.getCurrentLabel(); 
						if(lab=='S'||lab=='B'){
//System.out.println("-> "+reader.getCurrentChromosome()+" "+reader.getCurrentStart()+" "+j);
							scPos.add(preds[l][j]);
						}else if(lab=='U'){
							scNeg.add(preds[l][j]);
						}
						
					}else{
						break;
					}
					j++;					
				}while(reader.readNextFeatureVector());
			}

			
			double prev = ToolBox.sum(weights[1])*0.15;
			double perc2 = 1.0 - prev/(double)scNeg.length();
			
			double th = ToolBox.percentile(scPos.toArray(), perc);
			double th2 = ToolBox.percentile(scNeg.toArray(), perc2);
			
			th = Math.max(th, th2);
			System.out.println("threshold: "+th);
			
			/*HashSet<Integer> li = new HashSet<>(); 

			int k=0;
			for(int l=0;l<itChroms.size();l++){
				String chr = itChroms.get(l);
//System.out.println(chr);
				//double[] preds = pred.predict(chr, sizes.get(chr));
				for(int j=0;j<preds.length;j++,k++){
					if(preds[l][j] >= th){
						li.add(k);
					}
				}
			}*/
//System.out.println("predicted");
			LinkedList<Sequence> seqs = new LinkedList<>();
			
			reader.reset();
			for(int l=0;l<itChroms.size();l++){
				String chr = itChroms.get(l);
				
				reader.findChr(chr);
				int j=0;
				int size = sizes.get(chr);
				do{
					if(j<size){
						
						if(preds[l][j] >= th && reader.getCurrentLabel()=='U'){
//System.out.println("-> "+reader.getCurrentChromosome()+" "+reader.getCurrentStart()+" "+k);
							Sequence seq = reader.getCurrentSequence();
							seqs.add(seq);
						}						
						
					}else{
						break;
					}
					j++;					
				}while(reader.readNextFeatureVector());
			}
			
			
			bg = FeatureReader.replaceNaN( DataSet.union(bg,new DataSet("",seqs)) );
			double[] temp = new double[weights[1].length+seqs.size()];
			Arrays.fill(temp, 1.0);
			System.arraycopy(weights[1], 0, temp, 0, weights[1].length);
			weights[1] = temp;
			
			cl = train(new DataSet[]{fg,bg},weights);
			clList.add(cl);
			
		}
		
		return clList.toArray(new GenDisMixClassifier[0]);
		
	}
	
	
	
}
