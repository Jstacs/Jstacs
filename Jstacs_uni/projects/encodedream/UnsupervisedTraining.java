package projects.encodedream;
import java.io.FileNotFoundException;
import java.io.IOException;
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
import java.util.Random;

import de.jstacs.algorithms.optimization.Optimizer;
import de.jstacs.classifiers.differentiableSequenceScoreBased.OptimizableFunction.KindOfParameter;
import de.jstacs.classifiers.differentiableSequenceScoreBased.gendismix.GenDisMixClassifier;
import de.jstacs.classifiers.differentiableSequenceScoreBased.gendismix.GenDisMixClassifierParameterSet;
import de.jstacs.classifiers.differentiableSequenceScoreBased.gendismix.LearningPrinciple;
import de.jstacs.classifiers.differentiableSequenceScoreBased.logPrior.DoesNothingLogPrior;
import de.jstacs.data.DataSet;
import de.jstacs.data.EmptyDataSetException;
import de.jstacs.data.WrongAlphabetException;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.sequenceScores.statisticalModels.differentiable.continuous.GaussianNetwork;
import de.jstacs.utils.DoubleList;
import de.jstacs.utils.Pair;
import de.jstacs.utils.ToolBox;
import projects.dimont.Interpolation;

public class UnsupervisedTraining {
	
	
	public enum Init{
		MOTIF,
		DNASE,
		BOTH
	};
	
	public enum Select{
		ALTERNATE,
		RANDOM,
		FULL
	}
	
	private FeatureReader reader;
	private int threads;
	private HashMap<String,Integer> sizes;
	private Init init;
	private Select select;
	
	public UnsupervisedTraining(FeatureReader reader, int threads, HashMap<String,Integer> sizes, Init init, Select select){
		this.reader = reader;
		this.threads = threads;
		this.sizes = sizes;
		this.init = init;
		this.select = select;
	}
	
	
	public GenDisMixClassifier train(DataSet[] data, double[][] weights) throws Exception{
		
		GaussianNetwork gn = new GaussianNetwork(new int[data[0].getElementLength()][0]);
		
		GenDisMixClassifierParameterSet params = new GenDisMixClassifierParameterSet(data[0].getAlphabetContainer(), gn.getLength(), Optimizer.QUASI_NEWTON_BFGS, 1E-6, 1E-6, 1E-4, false, KindOfParameter.PLUGIN, true, threads);
		
		GenDisMixClassifier cl = new GenDisMixClassifier(params, DoesNothingLogPrior.defaultInstance, LearningPrinciple.MCL, gn,gn);
		
		cl.train(data, weights);
		
		return cl;
		
	}
	
	
	
	public GenDisMixClassifier[] iterativeTraining(int iterations, LinkedList<String> trainChroms, double frac, double cons) throws Exception{
		
		
		LinkedList<Sequence> full = new LinkedList<>();
		LinkedList<Sequence> motifs = new LinkedList<>();
		LinkedList<Sequence> dnase = new LinkedList<>();
		
		Pair<double[], double[][]> pair = getInitialWeights(trainChroms,frac,full,dnase,motifs);
		double[] vals = pair.getFirstElement();
		double[][] weights = pair.getSecondElement(); 
		
		DataSet fullDS = FeatureReader.replaceNaN( new DataSet("", full) );
		DataSet dnaseDS = FeatureReader.replaceNaN( new DataSet("", dnase) );
		DataSet motifsDS = FeatureReader.replaceNaN( new DataSet("", motifs) );
		
		
		LinkedList<GenDisMixClassifier> clList = new LinkedList<>();
		
		GenDisMixClassifier cl = train(new DataSet[]{fullDS,fullDS}, weights);
		clList.add(cl);
		
		DataSet curr = null;
		curr = getCurr(fullDS,dnaseDS,motifsDS,curr);
		
		GenDisMixClassifier currCl = cl;
		if(select != Select.FULL){
			currCl = train(new DataSet[]{curr,curr}, weights);
		}
		
		
		for(int i=0;i<iterations;i++){
			
			vals = updateVals(vals,currCl,curr,cons);
			weights = updateWeights(vals,frac);
			
			cl = train(new DataSet[]{fullDS,fullDS}, weights);
			clList.add(cl);
			
			curr = getCurr(fullDS,dnaseDS,motifsDS,curr);
			
			currCl = cl;
			if(select != Select.FULL){
				currCl = train(new DataSet[]{curr,curr}, weights);
			}
			
		}		
		
		return clList.toArray(new GenDisMixClassifier[0]);
		
	}

	
	
	private DataSet getCurr(DataSet fullDS, DataSet dnaseDS, DataSet motifsDS, DataSet curr) throws EmptyDataSetException, WrongAlphabetException {
		
		Random r = new Random(127);
		
		if(select == Select.ALTERNATE){
			if(curr == null){
				if(init == Init.MOTIF || init == Init.BOTH){
					return dnaseDS;
				}else{
					return motifsDS;
				}
			}else if(curr == dnaseDS){
				return motifsDS;
			}else{
				return dnaseDS;
			}
		}else if(select == Select.FULL){
			return fullDS;
		}else{
			
			int[] starts = new int[fullDS.getElementLength()];
			for(int i=0;i<starts.length;i++){
				starts[i] = i;
			}
			
			int[] starts2 = new int[starts.length/2];
			for(int i=0;i<starts.length/2;i++){
				int idx = r.nextInt(starts.length);
				int temp = starts[i];
				starts[i] = starts[idx];
				starts2[i] = starts[idx];
				starts[idx] = starts[i];
			}
			int[] lengths = new int[starts2.length];
			Arrays.fill(lengths, 1);
			
			Sequence[] seqs = new Sequence[fullDS.getNumberOfElements()];
			for(int i=0;i<seqs.length;i++){
				Sequence seq = fullDS.getElementAt(i);
				seqs[i] = seq.getCompositeSequence(starts2, lengths);
			}
			
			return new DataSet("",seqs);
			
		}
	}


	private double[][] updateWeights(double[] vals, double frac) throws Exception {
		double[] weights = Interpolation.getWeight(null, vals, frac, Interpolation.RANK_LOG);
		return new double[][]{weights,Interpolation.getBgWeight(weights)};
	}


	private static double[] updateVals(double[] vals, GenDisMixClassifier cl, DataSet curr, double cons) throws Exception {
		double[] temp = ToolBox.zscore(vals);
		
		double[] sc = cl.getScores(curr);
		sc = ToolBox.zscore(sc);
		
		for(int i=0;i<sc.length;i++){
				sc[i] = sc[i] + temp[i]*cons;
		}

		return sc;
	}
	
	

	private Pair<double[], double[][]> getInitialWeights(LinkedList<String> trainChroms, double frac, LinkedList<Sequence> full, LinkedList<Sequence> dnase, LinkedList<Sequence> motifs) throws Exception {
		
		reader.reset();
		
		DoubleList motifScores = new DoubleList();
		DoubleList dnaseScores = new DoubleList();
		
		for(int l=0;l<trainChroms.size();l++){
			String chr = trainChroms.get(l);
			
			reader.findChr(chr);
			int j=0;
			int size = sizes.get(chr);
			do{
				if(j<size){
					
					if(init == Init.BOTH || init == Init.MOTIF){
						motifScores.add(reader.getCurrentMotifMax(0));
					}
					
					if(init == Init.BOTH || init == Init.DNASE){
						dnaseScores.add(reader.getCurrentDNaseMin());
					}
					
					full.add(reader.getCurrentSequence());
					dnase.add(reader.getCurrentDNaseSequence());
					motifs.add(reader.getCurrentMotifsSequence());
					
				}else{
					break;
				}
				j++;					
			}while(reader.readNextFeatureVector());
		}
		
		
		double[] ms = motifScores.toArray();
		double[] ds = dnaseScores.toArray();
		
		double[] vals = null;
		
		if(init == Init.MOTIF){
			vals = ms;
		}else if(init == Init.DNASE){
			vals = ds;
		}else{
			
			ToolBox.zscore(ms);
			ToolBox.zscore(ds);
			double mi = ToolBox.min(ms);
			for(int i=0;i<ms.length;i++){
				ms[i] -= mi;
			}
			mi = ToolBox.min(ds);
			for(int i=0;i<ds.length;i++){
				ds[i] -= mi;
			}
			
			vals = ms;
			for(int i=0;i<vals.length;i++){
				vals[i] += ds[i];
			}
		}
		
		
		
		return new Pair<double[],double[][]>(vals,updateWeights(vals, frac));
		
	}
	
	
	
}
