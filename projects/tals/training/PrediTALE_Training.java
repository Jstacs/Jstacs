package projects.tals.training;

import de.jstacs.algorithms.optimization.ConstantStartDistance;
import de.jstacs.algorithms.optimization.Optimizer;
import de.jstacs.algorithms.optimization.termination.SmallDifferenceOfFunctionEvaluationsCondition;
import de.jstacs.classifiers.differentiableSequenceScoreBased.OptimizableFunction.KindOfParameter;
import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.DNADataSet;
import de.jstacs.data.DataSet;
import de.jstacs.data.alphabets.DiscreteAlphabet;
import de.jstacs.io.FileManager;
import de.jstacs.utils.SafeOutputStream;
import projects.tals.RVDSequence;
import projects.tals.linear.LF0Conditional;
import projects.tals.linear.LFModularConditional9C;
import projects.tals.linear.LFPosition_mixture;
import projects.tals.linear.LFSpecificity_parallel_cond9C;
import projects.tals.linear.MSDFunction;

public class PrediTALE_Training {

	public static void main(String[] args) throws Exception {
		String TrainData=args[0];
		int numberOfStarts=50;
	
		
		AlphabetContainer alphabet12 = new AlphabetContainer(new DiscreteAlphabet(false, "A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V"));
		AlphabetContainer alphabet13 = new AlphabetContainer(new DiscreteAlphabet(false, "A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V","*"));
		
		String[] sepRVDs = new String[]{"HD", "NN", "NG", "NI"};
		
		double con = 0.1;
		
		String[] firstPos = new String[]{"D","N","G","I"};
		
		LFSpecificity_parallel_cond9C spec = new LFSpecificity_parallel_cond9C(alphabet12, alphabet13, RVDSequence.getContainerRVD(alphabet12, alphabet13), firstPos, sepRVDs, new String[]{"HD","NN","NG","HG","NI","NK"});
		LFPosition_mixture pos = new LFPosition_mixture();

		String FileXML = "model.xml";
		
		DataSet ds = new DNADataSet(TrainData,'>', new CachingRVDReferenceSequenceAnnotationParser("seq", ":", ",", alphabet12, alphabet13));
		
		int threads = 8;
		
		DataSet[] dss = MSDFunction.splitByTagAndSort(threads, ds, "group", "signal", "globalWeight",false, true);
		double[][] gw = new double[dss.length][];
		for(int i=0;i<dss.length;i++){
			gw[i] = new double[dss[i].getNumberOfElements()];
			for(int j=0;j<gw[i].length;j++){
				gw[i][j] = Double.parseDouble(dss[i].getElementAt(j).getSequenceAnnotationByType("globalWeight", 0).getIdentifier());
			}
		}
		
		System.out.println("groups: "+dss.length);
		
		
				LFModularConditional9C model = new LFModularConditional9C(
						new LF0Conditional(RVDSequence.getContainerRVD(alphabet12, alphabet13), new String[]{"HD","NN","NG","NI","NS"},new double[]{0.05,0.3,0.05,0.6}),
						spec,
						pos, 
						dss.length);
				
							
		double bestResult=Double.POSITIVE_INFINITY;
		double[] bestVals=null;
		MSDFunction fun =null;
		for(int i=0;i<numberOfStarts;i++){
			System.out.println("Start: "+(i+1));
			model.initializeFunctionRandomly(false);
			
			System.out.println(model);
			
			fun = new MSDFunction(con, false, model, threads, dss, gw, dss.length*2, "signal");
			
			fun.reset();
			
			double[] vals = fun.getParameters(KindOfParameter.LAST);
			
			Optimizer.optimize(Optimizer.QUASI_NEWTON_BFGS, fun, vals, new SmallDifferenceOfFunctionEvaluationsCondition(1E-10), 1E-10, new ConstantStartDistance(1E-10), SafeOutputStream.getSafeOutputStream(System.out));
			
			fun.setParams(vals);
			
			double aktResult=fun.evaluateFunction(vals);
			if(aktResult<bestResult){
				bestResult=aktResult;
				bestVals=vals;
			}
			System.out.println("current result (Start: "+(i+1)+"):"+ aktResult);
			System.out.println(model);
		}
	  fun.setParams(bestVals);
	  
	  System.out.println(fun.evaluateFunction(bestVals));
		
	  System.out.println(model);
		
		
	  FileManager.writeFile(FileXML, model.toXML());
	  
}


	

}
