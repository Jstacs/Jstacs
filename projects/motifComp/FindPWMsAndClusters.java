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

package projects.motifComp;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.StringReader;
import java.text.DecimalFormat;
import java.util.AbstractMap.SimpleEntry;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Date;
import java.util.HashMap;
import java.util.LinkedList;

import de.jstacs.DataType;
import de.jstacs.clustering.distances.DeBruijnMotifComparison;
import de.jstacs.clustering.hierachical.ClusterTree;
import de.jstacs.data.DataSet;
import de.jstacs.data.alphabets.DNAAlphabetContainer;
import de.jstacs.io.FileManager;
import de.jstacs.io.NonParsableException;
import de.jstacs.io.SparseStringExtractor;
import de.jstacs.io.XMLParser;
import de.jstacs.parameters.FileParameter;
import de.jstacs.parameters.ParameterSet;
import de.jstacs.parameters.SelectionParameter;
import de.jstacs.parameters.SimpleParameter;
import de.jstacs.parameters.SimpleParameterSet;
import de.jstacs.parameters.validation.NumberValidator;
import de.jstacs.results.CategoricalResult;
import de.jstacs.results.ListResult;
import de.jstacs.results.PlotGeneratorResult;
import de.jstacs.results.Result;
import de.jstacs.results.ResultSet;
import de.jstacs.results.ResultSetResult;
import de.jstacs.sequenceScores.statisticalModels.StatisticalModel;
import de.jstacs.sequenceScores.statisticalModels.trainable.PFMWrapperTrainSM;
import de.jstacs.tools.JstacsTool;
import de.jstacs.tools.ProgressUpdater;
import de.jstacs.tools.Protocol;
import de.jstacs.tools.ToolParameterSet;
import de.jstacs.tools.ToolResult;
import de.jstacs.tools.JstacsTool.ResultEntry;
import de.jstacs.tools.ui.galaxy.Galaxy;
import de.jstacs.tools.ui.galaxy.MultilineSimpleParameter;
import de.jstacs.utils.ComparableElement;
import de.jstacs.utils.PFMComparator;
import de.jstacs.utils.Pair;
import de.jstacs.utils.SeqLogoPlotter.SeqLogoPlotGenerator;
import de.jstacs.utils.ToolBox;

public class FindPWMsAndClusters implements JstacsTool{

	private static int n = 8;
	
	public static void main(String[] args) throws Exception {
		
		Galaxy gal = new Galaxy(" -Xms512M -Xmx2G", false, new FindPWMsAndClusters());
		
		gal.run(args);
		
	}

	private ClusterTree<StatisticalModel>[] trees;
	private HashMap<String,String[]> expMap;
	
	
	public FindPWMsAndClusters() throws IOException, NonParsableException{
		StringBuffer buffer = FileManager.readInputStream(FindPWMsAndClusters.class.getClassLoader().getResourceAsStream( "projects/motifComp/data/clusters.xml" ) );
		trees = (ClusterTree<StatisticalModel>[]) XMLParser.extractObjectForTags(buffer, "trees");
		
		buffer = FileManager.readInputStream(FindPWMsAndClusters.class.getClassLoader().getResourceAsStream( "projects/motifComp/data/encode_ids.txt" ) );
		
		expMap = new HashMap<String, String[]>();
		
		String[] lines = buffer.toString().split("\n");
		for(int i=1;i<lines.length;i++){
			String[] parts = lines[i].split("\t");
			expMap.put(parts[0], parts);
		}
		
	}
	
	private StatisticalModel getModel(DataSet ds) throws CloneNotSupportedException{
		double[][] pfm = PFMComparator.getPFM(ds);
		return new PFMWrapperTrainSM(DNAAlphabetContainer.SINGLETON, "motif estimated from data set", pfm, 4.0);
	}
	
	private ResultSet find(double t, Pair<String, double[]>[] profiles, LinkedList<double[][]> motifs) throws Exception {
		LinkedList<Result> ress = new LinkedList<Result>();
		for(int i=0;i<profiles.length;i++){
			LinkedList<Result> res = find(t, profiles[i].getFirstElement(),profiles[i].getSecondElement());
			if(motifs.get(i) != null){
				res.addFirst(new PlotGeneratorResult("Motif", "Sequence logo of "+profiles[i].getFirstElement(), 
						new SeqLogoPlotGenerator(motifs.get(i), 200), true));
			}
			ress.add(new ResultSetResult("Matches for "+profiles[i].getFirstElement(), "", null, new ResultSet(res.toArray(new Result[0]))));
		}
		return new ResultSet(ress.toArray(new Result[0]));
	}

	private LinkedList<Result> find(double t, String name, double[] profile) throws Exception {
		
		double l = profile.length;
		int n = (int) Math.round( Math.log(l)/Math.log(4.0) );
		if(Math.pow(4.0, n) != l){
			throw new RuntimeException();
		}
		
		LinkedList<ComparableElement<MotifMatch,Double>> matches = new LinkedList<ComparableElement<MotifMatch,Double>>();
		for(int i=0;i<trees.length;i++){
			StatisticalModel[] mods = trees[i].getClusterElements();
			double mean = 0.0;
			LinkedList<ComparableElement<String, Double>> locMatch = new LinkedList<ComparableElement<String,Double>>();
			for(int j=0;j<mods.length;j++){
				double[][] prof1 = DeBruijnMotifComparison.getProfilesForMotif( mods[j], n, false, false );
				Pair<Integer,Double> fwd = DeBruijnMotifComparison.compare( prof1[0], profile, mods[j].getLength() );
				//	System.out.println("rev:");
				double[][] prof2 = DeBruijnMotifComparison.getProfilesForMotif( mods[j], n, true, false );
				Pair<Integer,Double> rev = DeBruijnMotifComparison.compare( prof2[0], profile, mods[j].getLength() );
				
				double v = Math.max(fwd.getSecondElement(), rev.getSecondElement());
				
				if(v > t){
					locMatch.add(new ComparableElement<String, Double>(mods[j].toString(), v));
				}
				
				mean += v;
			}
			mean /= mods.length;
			
			if(locMatch.size() > 0){
				matches.add(new ComparableElement<FindPWMsAndClusters.MotifMatch, Double>(new MotifMatch(trees[i], locMatch.toArray(new ComparableElement[0])), mean));
			}
			
		}
		
		ComparableElement<MotifMatch,Double>[] m2 = matches.toArray(new ComparableElement[0]);
		Arrays.sort(m2);
		
		LinkedList<Result> ress = new LinkedList<Result>();
		
		for(int i=m2.length-1;i>=0;i--){
			
			PlotGeneratorResult pgr = new PlotGeneratorResult("Cluster tree", "", new MotifTreePlotter.MotifTreePlotGenerator(m2[i].getElement().tree, 100, n) , true);
			
			LinkedList<ResultSet> sets = new LinkedList<ResultSet>();
			
			for(int j=m2[i].getElement().matches.length-1;j>=0;j--){
				
				String id = m2[i].getElement().matches[j].getElement();
				double val = m2[i].getElement().matches[j].getWeight();
				String valStr = DecimalFormat.getInstance().format(val);
				
				String key = id.substring(id.indexOf("(")+1);
				key = key.substring(0, key.indexOf("-"));
				
				String[] enc = expMap.get(key);
				
				ResultSet res = new ResultSet(new Result[]{
						new CategoricalResult("Motif", "", id),
						new CategoricalResult("Similarity", "", valStr),
						//new CategoricalResult("ENCODE ID", "", enc[0]),
						new CategoricalResult("Type", "", enc[1]),
						new CategoricalResult("Target","", enc[2]),
						new CategoricalResult("Description", "", enc[3]),
						new CategoricalResult("Lab", "", enc[4]),
						new CategoricalResult("Link", "", "<a href=\""+enc[5]+"\" target=\"_blank\">"+enc[5]+"</a>"),
				});
				sets.add(res);
				
			}
			
			ListResult list = new ListResult("List of motif matches", "", null, sets.toArray(new ResultSet[0]));
			
			ResultSetResult rsr = new ResultSetResult("Matching cluster "+(m2.length-i), "", null,
					new ResultSet(new Result[]{pgr,list}));
			
			ress.add(rsr);
			
		}
		
		return ress;
	}

	private static class MotifMatch{
		
		private ClusterTree<StatisticalModel> tree;
		private ComparableElement<String,Double>[] matches;
		
		public MotifMatch(ClusterTree<StatisticalModel> tree, ComparableElement<String,Double>[] matches){
			this.tree = tree;
			this.matches = matches.clone();
			Arrays.sort(this.matches);
		}
		
	}

	@Override
	public ToolParameterSet getToolParameters() {
		try{
		MultilineSimpleParameter pwms = new MultilineSimpleParameter("PWMs/PFMs", "PWMs/PFMs in Jaspar format. PWMs contain nucleotide probabilities and PFMs contain nucleotide frequencies.", true);
		
		FileParameter data = new FileParameter("Data set", "A data set of aligned binding sites", "fasta", true);
		FileParameter profiles = new FileParameter("Score profiles", "Score profiles in pseudo-FastA format", "txt", true);
		
		SelectionParameter sp = new SelectionParameter(DataType.PARAMETERSET, new String[]{"PFMs/PWMs","Binding sites","Score profiles"}, 
				new ParameterSet[]{
				new SimpleParameterSet(pwms),
				new SimpleParameterSet(data),
				new SimpleParameterSet(profiles)
		}, 
				"Motif source", "You may specify the query motifs as PWMs/PFMs/PSSMs, "
						+ "by a data set of aligned binding sites that are used to build a PWM, "
						+ "or by a score profile computed on one of the supplied de Bruijn sequences.", true);
		
		SimpleParameter t = new SimpleParameter(DataType.DOUBLE, "Similarity threshold", "The threshold on the correlation of score profiles, between 0.5 and 1.0.", true, new NumberValidator<Double>(0.5, 1.0), 0.9);
		
		return new ToolParameterSet(getShortName(),sp,t);
		
		}catch(Exception e){
			e.printStackTrace();
			throw new RuntimeException();
		}
	}

	@Override
	public ToolResult run(ToolParameterSet parameters, Protocol protocol,
			ProgressUpdater progress, int threads) throws Exception {
		
		LinkedList<double[][]> motifs = new LinkedList<double[][]>();
		
		Pair<String,double[]>[] profiles = getProfiles(parameters,protocol, motifs);
		
		double t = (Double) parameters.getParameterAt(1).getValue();
		
		protocol.append("Searching for matches of "+profiles.length+" motifs.\n");
		
		ResultSet all = find(t, profiles, motifs);
		
		protocol.append("Finished.\n");
		
		return new ToolResult("Motif matches", "", null, all, parameters, getToolName(), new Date());
	}

	private Pair<String, double[]>[] getProfiles(ParameterSet parameters,
			Protocol protocol, LinkedList<double[][]> motifs) throws Exception {
		SelectionParameter sp = (SelectionParameter) parameters.getParameterAt(0);
		if(sp.getSelected() == 2){
			String content = ((FileParameter)((ParameterSet)sp.getValue()).getParameterAt(0)).getFileContents().getContent();
			String[] parts = content.split("\n");
			String id = null;
			LinkedList<Pair<String,double[]>> profs = new LinkedList<Pair<String,double[]>>();
			for(int i=0;i<parts.length;i++){
				if(parts[i].startsWith(">")){
					id = parts[i].substring(1).trim();
				}else{
					String[] subs = parts[i].split("\\s");
					double[] vals = new double[subs.length];
					for(int j=0;j<vals.length;j++){
						vals[j] = Double.parseDouble(subs[j]);
					}
					profs.add(new Pair<String, double[]>(id, vals));
				}
			}
			motifs.add(null);
			protocol.append("Loaded "+profs.size()+" score profiles.\n");
			return profs.toArray(new Pair[0]);
		}else if(sp.getSelected() == 0 || sp.getSelected() == 1){//PWM
			
			ArrayList<SimpleEntry<String, double[][]>> pwms = null;

			if(sp.getSelected() == 1){
				pwms = new ArrayList<SimpleEntry<String,double[][]>>();
				String content = ((FileParameter)((ParameterSet)sp.getValue()).getParameterAt(0)).getFileContents().getContent();
				StringReader reader = new StringReader(content);
				
				DataSet data = new DataSet(DNAAlphabetContainer.SINGLETON, new SparseStringExtractor(reader, '>', "", null));
				
				double[][] pfm = PFMComparator.getPFM(data);
				pwms.add(new SimpleEntry<String, double[][]>("PFM build from data",pfm));
				protocol.append("Built PWM from "+data.getNumberOfElements()+" sequences of length "+data.getElementLength()+".\n");
			}else{

				String content = (String) ((ParameterSet)sp.getValue()).getParameterAt(0).getValue();
			//	System.out.println("\ncontent: >>"+content+"<<");
				pwms = PFMComparator.readPFMsFromJasparFastA( new BufferedReader( new StringReader(content) ) );
				protocol.append("Loaded "+pwms.size()+" PWMs/PFM from Jaspar format.\n");
			}
			
			Pair<String,double[]>[] profiles = new Pair[pwms.size()];
			for(int i=0;i<profiles.length;i++){
				String name = pwms.get(i).getKey();
				double[][] pwm = pwms.get(i).getValue();
				
				boolean isPWM = true;
				for(int j=0;j<pwm.length;j++){
					if( ToolBox.sum(pwm[j]) > 1.0+1E-6 ){
						isPWM = false;
						break;
					}
				}
				PFMWrapperTrainSM model =  new PFMWrapperTrainSM(DNAAlphabetContainer.SINGLETON, name, pwm, isPWM ? 0.0 : 4.0);
				
				motifs.add(model.getPWM());
				
				double[] prof = DeBruijnMotifComparison.getProfilesForMotif(model, n, false, false)[0];
				protocol.append("Determined score profile for "+name+".\n");
				/*PrintWriter wr = new PrintWriter("/Users/dev/Downloads/scoretest.txt");
				wr.println(">test");
				for(int j=0;j<prof.length;j++){
					if(j>0){
						wr.print(" ");
					}
					wr.print(prof[j]);
				}
				wr.println();
				wr.close();*/
				
				profiles[i] = new Pair<String, double[]>(name, prof); 
			}
			return profiles;
		}else{
			return null;
		}
	}

	@Override
	public String getToolName() {
		return "DBcorrDB";
	}

	@Override
	public String getShortName() {
		return "dbcorrdb";
	}

	@Override
	public String getDescription() {
		return "search a data base of motifs by similarity of score profiles on de Bruijn sequences";
	}

	@Override
	public String getHelpText() {
		try {
			return FileManager.readInputStream( FindPWMsAndClusters.class.getClassLoader().getResourceAsStream( "projects/motifComp/FindPWMsAndClusters.txt" ) ).toString();
		} catch ( IOException e ) {
			e.printStackTrace();
			return "";
		}
	}
	
	@Override
	public String getToolVersion() {
		return "1.0";
	}

	@Override
	public ResultEntry[] getDefaultResultInfos() {
		return null;
	}

	@Override
	public ToolResult[] getTestCases(String path) {
		return null;
	}

	@Override
	public void clear() {		
	}

	@Override
	public String[] getReferences() {
		return null;
	}
}
