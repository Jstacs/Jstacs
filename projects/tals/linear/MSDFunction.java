package projects.tals.linear;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;

import de.jstacs.algorithms.optimization.DimensionException;
import de.jstacs.algorithms.optimization.EvaluationException;
import de.jstacs.classifiers.differentiableSequenceScoreBased.AbstractMultiThreadedOptimizableFunction;
import de.jstacs.classifiers.differentiableSequenceScoreBased.OptimizableFunction.KindOfParameter;
import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.DataSet;
import de.jstacs.data.EmptyDataSetException;
import de.jstacs.data.WrongAlphabetException;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.data.sequences.annotation.ReferenceSequenceAnnotation;
import de.jstacs.data.sequences.annotation.SequenceAnnotation;
import de.jstacs.sequenceScores.differentiable.DifferentiableSequenceScore;
import de.jstacs.utils.ComparableElement;
import de.jstacs.utils.DoubleList;
import de.jstacs.utils.IntList;
import de.jstacs.utils.ToolBox;

public class MSDFunction extends AbstractMultiThreadedOptimizableFunction {


	private DifferentiableSequenceScore[] scores;

	private double[][][] yi;
	
	private double[][] grads;
	
	private double[] vals;
	
	private IntList[] indices;
	private DoubleList[] partDers;
	
	private double con;
	private boolean laplace;
	private int penaltyOff;
	
	private String sortTag;
	
	public MSDFunction(double con, boolean laplace, DifferentiableSequenceScore score, int threads, DataSet[] data, double[][] weights, int penaltyOff, String sortTag) throws IllegalArgumentException {
		super(threads, data, weights, false, false);
		this.scores = new DifferentiableSequenceScore[threads];
		this.scores[0] = score;
		precomputeIndexes(sortTag);
		//numEvals = 0;
		this.con = con;
		this.laplace = laplace;
		this.penaltyOff = penaltyOff;
		this.sortTag = sortTag;
	}

	@Override
	public int getDimensionOfScope() {
		return this.scores[0].getNumberOfParameters();
	}

	@Override
	protected void evaluateGradientOfFunction(int index, int startClass, int startSeq, int endClass, int endSeq) {
		
		if(startSeq != 0 || endSeq != data[endClass].getNumberOfElements()){
			throw new RuntimeException();
		}
		
		Arrays.fill(grads[index], 0);
		
		for(int cl = startClass; cl<=endClass; cl++){
			int start=0, end = data[cl].getNumberOfElements();
			
			for(int i=start;i<end;i++){
				indices[index].clear();
				partDers[index].clear();
				
				
				double s = scores[index].getLogScoreAndPartialDerivation(data[cl].getElementAt(i), indices[index], partDers[index]);
				
				double v = 2*(s-yi[index][cl-startClass][i])*weights[cl][i];
				
				for(int j=0;j<indices[index].length();j++){
					grads[index][ indices[index].get(j) ] += v * partDers[index].get(j);
				}
				
				
			}
		}
		
		
	}

	@Override
	protected double[] joinGradients() throws EvaluationException {
		for(int i=1;i<grads.length;i++){
			for(int j=0;j<grads[0].length;j++){
				grads[0][j] += grads[i][j];
			}
		}
		
		
		for(int i=0;i<params.length;i++){
			if(i>=penaltyOff){
				if(laplace){
					grads[0][i] += (params[i] < 0 ? -1 : 1)*con;
				}else{
					grads[0][i] += 2*params[i]*con;
				}
			}
			
		}
		
		for(int i=0;i<params.length;i++){
			grads[0][i] /= sum[cl];
			
		}
		
		
		//System.out.println("grad: "+Arrays.toString(grads[0]));
		return grads[0].clone();
	}

	@Override
	protected void evaluateFunction(int index, int startClass, int startSeq, int endClass, int endSeq)
			throws EvaluationException {
		
		if(startSeq != 0 || endSeq != data[endClass].getNumberOfElements()){
			throw new RuntimeException();
		}
		
		double val = 0;
		
	/*	double meanA = 0.0;
		double meanB = 0.0;
		double num = 0.0;*/
		
		for(int cl = startClass; cl<=endClass; cl++){
			int start=0, end=data[cl].getNumberOfElements();
			
			/*double meanScore = 0.0;
			double meanSq = 0.0;
			double meanY = 0.0;
			double meanSY = 0.0;
			double n = 0.0;*/
			
			for(int i=start;i<end;i++){
				double s = scores[index].getLogScoreFor( data[cl].getElementAt(i) );
				
				/*meanScore += s*weights[cl][i];
				meanY += yi[index][cl-startClass][i]*weights[cl][i];
				meanSq += s*s*weights[cl][i];
				meanSY += s*yi[index][cl-startClass][i]*weights[cl][i];
				n += weights[cl][i];*/
				
				double v = (s-yi[index][cl-startClass][i]);
				val += v*v*weights[cl][i];
				
			}
			

			/*meanScore /= n;
			meanY /= n;
			meanSq /= n;
			meanSY /= n;
			double myA = ( meanSY - meanScore*meanY ) / ( meanSq - meanScore*meanScore );
			double myB = meanY - myA*meanScore;
			//System.out.println(cl+" "+myA+" "+myB);
			if(!Double.isNaN(myA)){
				meanA += myA;
				meanB += myB;
				num ++;
			}*/

		}
		
		//System.out.println(index+" "+(meanA/num)+" "+(meanB/num));
		
		
		vals[index] = val;
		
	}

	@Override
	protected double joinFunction() throws EvaluationException, DimensionException {
		
		
		double val = ToolBox.sum(vals);
		
		for(int i=0;i<params.length;i++){
			if(i>=penaltyOff){
				if(laplace){
					val += Math.abs(params[i])*con;
				}else{
					val += params[i]*params[i]*con;
				}
			}
		}
		
		val /= sum[cl];
		
		//System.out.println(val);
		return val;
		
	}
	
	//assumed to be symmetric, i.e. w_i,j=w_j,i
	/*private double computeWeight(int cl, int i, int j, int numberOfElements){
		double w1 = Math.abs(i-j)/(double)numberOfElements;
		//double w2 = 1.0 - Math.min(i, j)/(double)numberOfElements;
		double w2 = Math.exp(-Math.min(i, j)/(double)numberOfElements*100.0);
		//double w2 = 1.0-Math.min(i, j)/(double)numberOfElements;
		//double w2 = 1.0 - Math.pow( Math.min(i, j)/(double)numberOfElements, 0.01);
		return w1*w2*weights[cl][i]*weights[cl][j];
	}*/
	
	/*private double getWeight2(int cl, int i, int j, int numberOfElements) {
		double w1 = Math.abs(i-j)/(double)numberOfElements;
		//double w2 = 1.0 - Math.min(i, j)/(double)numberOfElements;
		double w2 = Math.exp(-Math.min(i, j)/(double)numberOfElements*100.0);
		//double w2 = 1.0 - Math.pow( Math.min(i, j)/(double)numberOfElements, 0.01);
		return w1*w2*weights[cl][i]*weights[cl][j];
		//return 1E-4;
	}*/

	@Override
	protected void setThreadIndependentParameters() throws DimensionException {
		
	}

	

	@Override
	public void setDataAndWeights(DataSet[] data, double[][] weights) throws IllegalArgumentException {
		super.setDataAndWeights(data, weights);
		precomputeIndexes(sortTag);
	}

	public void precomputeIndexes(String sortTag) {
		if(worker != null){
			
			vals = new double[worker.length];
			indices = new IntList[worker.length];
			partDers = new DoubleList[worker.length];
			grads = new double[worker.length][getDimensionOfScope()];
			yi = new double[worker.length][][];
			
			for(int i=0;i<worker.length;i++){
				
				
				indices[i] = new IntList();
				partDers[i] = new DoubleList();
				
				
				int[] temp = worker[i].getIndices();
				int startClass = temp[0];
				int startSeq = temp[1];
				int endClass = temp[2];
				int endSeq = temp[3];

				if(startSeq != 0 || endSeq != data[endClass].getNumberOfElements()){
					throw new RuntimeException();
				}
				
				yi[i] = new double[endClass-startClass+1][];
				
				for(int j=startClass;j<=endClass;j++){
					yi[i][j-startClass] = new double[data[j].getNumberOfElements()];
					
					for(int k=0;k<data[j].getNumberOfElements();k++){
						Sequence seq = data[j].getElementAt(k);
						double y = Double.parseDouble(seq.getSequenceAnnotationByType(sortTag, 0).getIdentifier());
						yi[i][j-startClass][k] = y;
					}
				}
				
				
			}
			
		}
		
		
		
	}

	@Override
	protected void setParams(int index) throws DimensionException {
		
		this.scores[index].setParameters(params, 0);

	}

	@Override
	public void getParameters(KindOfParameter kind, double[] erg) throws Exception {
		double[] temp = this.scores[0].getCurrentParameterValues();
		System.arraycopy(temp, 0, erg, 0, temp.length);
	}

	@Override
	public void reset() throws Exception {
		for(int i=1;i<scores.length;i++){
			this.scores[i] = (DifferentiableSequenceScore) scores[0].clone();
		}
	}
	
	
	
	protected void prepareThreads() {
		double[] sizes = new double[data.length];

		
		for(int i=0;i<data.length;i++){
			sizes[i] = data[i].getNumberOfElements();
			sizes[i] = sizes[i]*sizes[i]*Math.sqrt(data[i].getAverageElementLength());
		}
		
		//System.out.println("sizes:"+Arrays.toString(sizes));
		
		double sum = ToolBox.sum(sizes);
		double part = sum/(double)worker.length;
		
		int startClass = 0;
		
		for(int i=0;i<worker.length;i++){
			
			double curr = 0;
			int endClass = startClass;
			curr = sizes[endClass];
			while( endClass < sizes.length-1 && curr+sizes[endClass+1] <= part ){
				curr += sizes[endClass+1];
				endClass++;
			}
			sum -= curr;
			//curr += sizes[endClass];
			
			if(i==worker.length-1){
				endClass = data.length-1;
			}
			
			/*System.out.println(i+":"+curr+" "+part);
			for(int j=startClass;j<=endClass;j++){
				System.out.println(" "+j+":"+data[j].getNumberOfElements()+" "+data[j].getAverageElementLength());
			}*/
			
			if( worker[i] != null ) {
				if( worker[i].isWaiting() ) {
					worker[i].setIndices(startClass, 0, endClass, data[endClass].getNumberOfElements());
				} else {
					stopThreads();
					throw new RuntimeException();
				}
			} else {
				worker[i] = new Worker( i, startClass, 0, endClass, data[endClass].getNumberOfElements() );
				worker[i].start();
			}
			
			startClass = endClass+1;
			
			part = sum/(double)(worker.length-i-1);
			
		}
		
		//System.out.println(Arrays.toString(worker));
	}
	
	
	private static double[] getMAD(double[] ws, double[] gw, boolean print){
		ws = ws.clone();
		gw = gw.clone();
		ToolBox.sortAlongWith(ws, gw);
		if(print) System.out.println(Arrays.toString(ws)+" "+Arrays.toString(gw));
		
		double fullSum = ToolBox.sum(gw);
		
		double median = 0.0;
		
		double sum = fullSum;
		double percentile = sum*0.5;
		for(int j=0;j<ws.length;j++){
			sum -= gw[j];
			if(print) System.out.println(sum+"<"+percentile+"?");
			if(sum*(1.0+1E-6) < percentile){
				double w1 = gw[j];
				double w2 = gw[j-1];
				median = (ws[j]*w1+ws[j-1]*w2)/(w1+w2);
				break;
			}
		}
		if(Double.isNaN(median) || Double.isInfinite(median)){
			median = 0.0;
		}
		double[] mads = new double[ws.length];
		for(int j=0;j<ws.length;j++){
			mads[j] = Math.abs(median - ws[j]);//-weight!
		}
		ToolBox.sortAlongWith(mads, gw);
		
		double mad = 1.0;
		
		sum = fullSum;
		for(int j=0;j<mads.length;j++){
			sum -= gw[j];
			if(sum*(1.0+1E-6) < percentile){
				double w1 = gw[j];
				double w2 = gw[j-1];
				mad = (mads[j]*w1 + mads[j-1]*w2)/(w1+w2);
				break;
			}
			
		}
		if(mad <= 1E-6 || Double.isNaN(mad) || Double.isInfinite(mad)){
			mad = 1.0;
		}
		
		return new double[]{median,mad};
		
	}
	
	public static DataSet[] splitByTagAndSort(int numThreads, DataSet data, String splitTag, String sortTag, String globalWeightTag, boolean filter, boolean normalize) throws EmptyDataSetException, WrongAlphabetException{
		
		HashMap<String, LinkedList<Sequence>> sets = new HashMap<String, LinkedList<Sequence>>();
		
		for(int i=0;i<data.getNumberOfElements();i++){
			Sequence seq = data.getElementAt(i);
			SequenceAnnotation sa = seq.getSequenceAnnotationByType(splitTag, 0);
			
			String key = "null";
			if(sa != null){
				key = seq.getSequenceAnnotationByType(splitTag, 0).getIdentifier();
			}
			if(!sets.containsKey(key)){
				sets.put(key, new LinkedList<Sequence>());
			}
			sets.get(key).add(seq);
		}
		
		DataSet[] ds = new DataSet[sets.keySet().size()]; 
		
		int i=0;
		Iterator<String> keys = sets.keySet().iterator();
		
		while(keys.hasNext()){
			String key = keys.next();
			
			Sequence[] seqs = sets.get(key).toArray(new Sequence[0]);
			
			ComparableElement<Sequence, Double>[] ws = new ComparableElement[seqs.length];
			double sum = 0.0;
			double sumsq = 0.0;
			double n = 0.0;
			double[] gws = new double[seqs.length];
			double[] lws = new double[seqs.length];
			for(int j=0;j<seqs.length;j++){
				double w = Double.parseDouble( seqs[j].getSequenceAnnotationByType(sortTag, 0).getIdentifier() );
				
				SequenceAnnotation an = seqs[j].getSequenceAnnotationByType(globalWeightTag, 0);
				
				double gw = Double.parseDouble( an.getIdentifier() );
				sum += w*gw;
				sumsq += w*w*gw;
				n += gw;
				ws[j] = new ComparableElement<Sequence, Double>(seqs[j], -w);
				
				gws[j] = gw;
				lws[j] = w;
			}
			
			
			sumsq /= n;
			double mean = sum/n;
			double sd = Math.sqrt(sumsq - mean*mean);
			
			if(n == 0){
				mean = 0.0;
				sd = 1.0;
			}
			if(sd <= 0){
				System.err.println(ws[0].getElement().getSequenceAnnotationByType(splitTag, 0).getIdentifier()+" "+mean+" "+sd);
				sd = 1.0;
			}
			
			Arrays.sort(ws);
			
			double min = -ws[ws.length-1].getWeight();
			double max = -ws[0].getWeight();
			
			/*double[] temp = getMAD(lws,gws,false);
			double median = temp[0];
			double mad = temp[1];
			if(sd/mad > 100){
				System.out.println(ws[0].getElement().getSequenceAnnotationByType(splitTag, 0).getIdentifier()+" "+min+" "+max+" ("+mean+" "+sd+") ("+median+" "+mad+")");
				getMAD(lws,gws,true);
			}
			
			double mmm = max-median;
			if(mmm < 1){
				System.out.println(ws[0].getElement().getSequenceAnnotationByType(splitTag, 0).getIdentifier()+" "+mmm+" "+max+" "+median);
			}*/
			
			
			for(int j=0;j<ws.length;j++){
				seqs[j] = ws[j].getElement();
				SequenceAnnotation mask = seqs[j].getSequenceAnnotationByType("mask", 0);
				
				int numPos = mask == null ? seqs[j].getLength() : (mask.getIdentifier().length() - mask.getIdentifier().replaceAll("X", "").length() );
				
				double lw = Double.parseDouble(seqs[j].getSequenceAnnotationByType(sortTag, 0).getIdentifier());
				
				double myMax = max;
				
				/*if(!normalize && seqs[j].getSequenceAnnotationByType(splitTag, 0).getIdentifier().startsWith("B")){
					lw -= min;
					myMax -= min;
				}*/
				
				if(normalize){
					lw = (lw - mean)/sd;
				}/*else{
					if(max == min){
						lw = 0;
					}else{
						lw = (max-lw)/(max-min)*Math.log(0.25)*numPos;
					}
				}*/
				//System.out.println(seqs[j].getSequenceAnnotationByType(splitTag, 0).getIdentifier()+" "+seqs[j]+" "+lw+" "+max+" "+min+" "+myMax);
				//lw /= sd;
				
				SequenceAnnotation mms = seqs[j].getSequenceAnnotationByType("mms", 0);
				
				if(mask!= null){
					seqs[j] = seqs[j].annotate(false, 
							new SequenceAnnotation("intgroup", i+""),
							(ReferenceSequenceAnnotation)seqs[j].getSequenceAnnotationByType("reference", 0),
							mask,
							seqs[j].getSequenceAnnotationByType(globalWeightTag, 0),
							new SequenceAnnotation( sortTag, lw+"" )
							);
				}else{
					seqs[j] = seqs[j].annotate(false, 
							new SequenceAnnotation("intgroup", i+""),
							(ReferenceSequenceAnnotation)seqs[j].getSequenceAnnotationByType("reference", 0),
							seqs[j].getSequenceAnnotationByType(globalWeightTag, 0),
							new SequenceAnnotation( sortTag, lw+"" )
							);
				}
				if(mms != null){
					seqs[j] = seqs[j].annotate(true, mms);
				}
			}
			
			
			
			
			if(filter){
				
				ArrayList<Sequence> list = new ArrayList<Sequence>();
				for(int j=0;j<seqs.length;j++){
					SequenceAnnotation mask = seqs[j].getSequenceAnnotationByType("mask", 0);
					if(mask == null || mask.getIdentifier().indexOf("X") < 0){
						ReferenceSequenceAnnotation an = (ReferenceSequenceAnnotation)seqs[j].getSequenceAnnotationByType("reference", 0);
						Sequence ref = an.getReferenceSequence();
						AlphabetContainer rvds = ref.getAlphabetContainer();
						int nmm = 0;
						for(int k=0;k<ref.getLength();k++){
							if(ref.discreteVal(k) == rvds.getCode(k, "HD")){
								if(seqs[j].discreteVal(k+1) != 1){
									nmm++;
								}
							}else if(ref.discreteVal(k) == rvds.getCode(k, "NI")){
								if(seqs[j].discreteVal(k+1) != 0){
									nmm++;
								}
							}else if(ref.discreteVal(k) == rvds.getCode(k, "NG")){
								if(seqs[j].discreteVal(k+1) != 3){
									nmm++;
								}
							}else if(ref.discreteVal(k) == rvds.getCode(k, "NN")){
								if(seqs[j].discreteVal(k+1) != 0 && seqs[j].discreteVal(k+1) != 2){
									nmm++;
								}
							}else{
								nmm = 0;//other RVD
								break;
							}
						}
						if(nmm <= 3){
							list.add(seqs[j]);
						}
						
					}else{
						list.add(seqs[j]);
					}
				}
				seqs = list.toArray(new Sequence[0]);
				
			}
			
			
			
			
			ds[i] = new DataSet("",seqs);
			
			i++;
		}
		
		if(numThreads > 1){

			double[] sizes = new double[ds.length];

			
			for(i=0;i<ds.length;i++){
				sizes[i] = ds[i].getNumberOfElements();
				sizes[i] = sizes[i]*sizes[i]*Math.sqrt(ds[i].getAverageElementLength());
			}
			
			
			
			
			int[] order = ToolBox.order(sizes, true);
			/*for(int j=0;j<order.length;j++){
				System.out.println(j+" "+sizes[order[j]]);
			}*/
			
			IntList[] lists = new IntList[numThreads];
			for(int j=0;j<lists.length;j++){
				lists[j] = new IntList();
			}
			double[] curr = new double[numThreads];
			
			for(int j=0;j<order.length;j++){
				double size = sizes[order[j]];
				int idx = ToolBox.getMinIndex(curr);
				lists[idx].add(order[j]);
				curr[idx] += size;
			}
			
			//System.out.println("curr:" + Arrays.toString(curr));
			
			DataSet[] ds2 = new DataSet[ds.length];
			
			for(int j=0,k=0;j<lists.length;j++){
				for(int l=0;l<lists[j].length();l++,k++){
					ds2[k] = ds[lists[j].get(l)];
				}
			}
			
			ds = ds2;
			
		}
		
		
		return ds;
		
	}
	
	
	

}
