package projects.tals.linear;

import java.text.NumberFormat;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;

import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.DataSet;
import de.jstacs.data.WrongAlphabetException;
import de.jstacs.data.alphabets.DNAAlphabetContainer;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.data.sequences.annotation.ReferenceSequenceAnnotation;
import de.jstacs.io.ArrayHandler;
import de.jstacs.io.NonParsableException;
import de.jstacs.io.XMLParser;
import de.jstacs.sequenceScores.differentiable.AbstractDifferentiableSequenceScore;
import de.jstacs.utils.DoubleList;
import de.jstacs.utils.IntList;
import projects.tals.RVDSequence;

public class LFSpecificity_parallel_cond9C extends AbstractDifferentiableSequenceScore {

	private int offSpecsThirteen, offSeparateSpecs, offConditionalSpecs;
	
	private int[] firstPosMap;
	private double[][] firstPosSpecs;
	
	private double[][] specsThirteen;
	
	private int[] conditionalMapRVD;
	private int[] conditionalMapTwelve;
	private double[][][] conditionalSpecs;
	
	private int[] separateMap;
	private double[][] separateSpecs;
	
	private AlphabetContainer twelve,thirteen,rvds;
	
	public LFSpecificity_parallel_cond9C( AlphabetContainer twelve,AlphabetContainer thirteen, AlphabetContainer rvds, String[] firstPosDifferent, String[] conditionalRVDs, String[] separateRVDs ) throws IllegalArgumentException, WrongAlphabetException {
		super(DNAAlphabetContainer.SINGLETON, 1);
		this.twelve = twelve;
		this.thirteen = thirteen;
		this.rvds = rvds;
		
		firstPosMap = new int[(int) thirteen.getAlphabetLengthAt(0)];
		Arrays.fill(firstPosMap, -1);
		firstPosSpecs = new double[firstPosDifferent.length][(int) this.alphabets.getAlphabetLengthAt(0)];
		int k=0;
		for(int i=0;i<firstPosDifferent.length;i++){
			firstPosMap[ (int) thirteen.getCode(0, firstPosDifferent[i]) ] = k;
			k++;
		}
		
		specsThirteen = new double[(int) thirteen.getAlphabetLengthAt(0)][(int) this.alphabets.getAlphabetLengthAt(0)];
		
		HashSet<String> condTwelves = new HashSet<>();
		for(int i=0;i<conditionalRVDs.length;i++){
			condTwelves.add(conditionalRVDs[i].substring(0, 1));
		}
		
		conditionalMapTwelve = new int[(int) twelve.getAlphabetLengthAt(0)];
		Arrays.fill(conditionalMapTwelve, -1);
		k=0;
		for(String condTwelve : condTwelves){
			conditionalMapTwelve[ (int) twelve.getCode(0, condTwelve) ] = k;
			k++;
		}
		
		conditionalMapRVD = new int[(int) rvds.getAlphabetLengthAt(0)];
		Arrays.fill(conditionalMapRVD, -1);
		k=0;
		for(String condRVD : conditionalRVDs){
			conditionalMapRVD[ (int) rvds.getCode(0, condRVD) ] = k;
			k++;
		}
		conditionalSpecs = new double[condTwelves.size()][conditionalRVDs.length][(int) this.alphabets.getAlphabetLengthAt(0)];
		
		separateMap = new int[(int) rvds.getAlphabetLengthAt(0)];
		Arrays.fill(separateMap, -1);
		k=0;
		for(String sepRVD : separateRVDs){
			separateMap[ (int) rvds.getCode(0, sepRVD) ] = k;
			k++;
		}
		separateSpecs = new double[separateRVDs.length][(int) this.alphabets.getAlphabetLengthAt(0)];
		
		offSpecsThirteen = firstPosSpecs.length*(firstPosSpecs.length > 0 ? firstPosSpecs[0].length : 0);
		offSeparateSpecs = offSpecsThirteen + specsThirteen.length * (specsThirteen.length > 0 ? specsThirteen[0].length : 0);
		offConditionalSpecs = offSeparateSpecs + separateSpecs.length * (separateSpecs.length > 0 ? separateSpecs[0].length : 0);
		
	}

	public LFSpecificity_parallel_cond9C(StringBuffer xml) throws NonParsableException {
		super(xml);
	}
	
	public LFSpecificity_parallel_cond9C clone() throws CloneNotSupportedException{
		LFSpecificity_parallel_cond9C clone = (LFSpecificity_parallel_cond9C) super.clone();
		
		clone.conditionalMapRVD = conditionalMapRVD.clone();
		clone.conditionalMapTwelve = conditionalMapTwelve.clone();
		clone.conditionalSpecs = ArrayHandler.clone(conditionalSpecs);
		clone.firstPosMap = firstPosMap.clone();
		clone.firstPosSpecs = ArrayHandler.clone(firstPosSpecs);
		clone.separateMap = separateMap.clone();
		clone.separateSpecs = ArrayHandler.clone(separateSpecs);
		clone.specsThirteen = ArrayHandler.clone(specsThirteen);
		
		return clone;
	}

	@Override
	public void initializeFunction(int index, boolean freeParams, DataSet[] data, double[][] weights) throws Exception {
		initializeFunctionRandomly(freeParams);
	}

	private void initVec(double[] vec){
		for(int i=0;i<vec.length;i++){
			vec[i] = 2.0*Math.random() - 1.0;
		}
	}
	
	private void initMatrix(double[][] mat){
		for(int i=0;i<mat.length;i++){
			initVec(mat[i]);
		}
	}
	
	@Override
	public void initializeFunctionRandomly(boolean freeParams) throws Exception {
		for(int i=0;i<conditionalSpecs.length;i++){
			initMatrix(conditionalSpecs[i]);
		}
		initMatrix(firstPosSpecs);
		initMatrix(separateSpecs);
		initMatrix(specsThirteen);
	}

	@Override
	public double getLogScoreAndPartialDerivation(Sequence seq, int start, IntList indices, DoubleList partialDer) {
		
		ReferenceSequenceAnnotation data_anno =(ReferenceSequenceAnnotation)seq.getSequenceAnnotationByType("reference", 0);
		RVDSequence rvd_seq=(RVDSequence)data_anno.getReferenceSequence();
		
		int index_rj = rvd_seq.discreteVal(start-1);//current RVD
		
		int index_dj=rvd_seq.discreteValThirteen( index_rj );//13 current
		int index_zj=rvd_seq.discreteValTwelve( index_rj );//12 current
		
		int index_rjp = -1;//previous RVD
		int index_djp=-1;//13 previous
		int index_zjp=-1;//12 previous
		if(start-2 >= 0){
			index_rjp = rvd_seq.discreteVal(start-2);
			index_djp=rvd_seq.discreteValThirteen( index_rjp );
			index_zjp=rvd_seq.discreteValTwelve( index_rjp );
		}
		
		
		
		double score = specsThirteen[index_dj][seq.discreteVal(start)];
		
		indices.add( offSpecsThirteen +  index_dj*specsThirteen[index_dj].length + seq.discreteVal(start));
		partialDer.add(1.0);
		
		if(separateMap[index_rj] >= 0){
			score += separateSpecs[separateMap[index_rj]][seq.discreteVal(start)];
			
			indices.add(offSeparateSpecs + separateMap[index_rj]*separateSpecs[separateMap[index_rj]].length + seq.discreteVal(start));
			partialDer.add(1.0);
		}
		
		if(start == 1 && firstPosMap[index_dj] >= 0){//first RVD and specific probability for that 13th AA
			
			score += firstPosSpecs[firstPosMap[index_dj]][seq.discreteVal(start)];
			
			indices.add(firstPosMap[index_dj]*firstPosSpecs[firstPosMap[index_dj]].length + seq.discreteVal(start));
			partialDer.add(1.0);
			
		}else if(index_zjp >= 0 && conditionalMapTwelve[index_zjp] >= 0 && conditionalMapRVD[index_rj] >= 0){
			
			score += conditionalSpecs[conditionalMapTwelve[index_zjp]][conditionalMapRVD[index_rj]][seq.discreteVal(start)];
			
			indices.add(offConditionalSpecs + 
					conditionalMapTwelve[index_zjp]*conditionalSpecs[conditionalMapTwelve[index_zjp]].length*conditionalSpecs[conditionalMapTwelve[index_zjp]][conditionalMapRVD[index_rj]].length +
					conditionalMapRVD[index_rj]*conditionalSpecs[conditionalMapTwelve[index_zjp]][conditionalMapRVD[index_rj]].length +
					seq.discreteVal(start));
			partialDer.add(1.0);
			
		}
		
		return score;
		
	}

	@Override
	public int getNumberOfParameters() {
		return offConditionalSpecs + conditionalSpecs.length * (conditionalSpecs.length > 0 ? conditionalSpecs[0].length * conditionalSpecs[0][0].length : 0);
	}

	@Override
	public double[] getCurrentParameterValues() throws Exception {
		double[] pars = new double[getNumberOfParameters()];
		int start = 0;
		for(int i=0;i<firstPosSpecs.length;i++){
			System.arraycopy(firstPosSpecs[i], 0, pars, start, firstPosSpecs[i].length);
			start += firstPosSpecs[i].length;
		}
		for(int i=0;i<specsThirteen.length;i++){
			System.arraycopy(specsThirteen[i], 0, pars, start, specsThirteen[i].length);
			start += specsThirteen[i].length;
		}
		for(int i=0;i<separateSpecs.length;i++){
			System.arraycopy(separateSpecs[i], 0, pars, start, separateSpecs[i].length);
			start += separateSpecs[i].length;
		}
		for(int i=0;i<conditionalSpecs.length;i++){
			for(int j=0;j<conditionalSpecs[i].length;j++){
				System.arraycopy(conditionalSpecs[i][j], 0, pars, start, conditionalSpecs[i][j].length);
				start += conditionalSpecs[i][j].length;
			}
		}
		return pars;
	}

	@Override
	public void setParameters(double[] params, int start) {
		for(int i=0;i<firstPosSpecs.length;i++){
			System.arraycopy(params, start, firstPosSpecs[i], 0, firstPosSpecs[i].length);
			start += firstPosSpecs[i].length;
		}
		for(int i=0;i<specsThirteen.length;i++){
			System.arraycopy(params, start, specsThirteen[i], 0, specsThirteen[i].length);
			start += specsThirteen[i].length;
		}
		for(int i=0;i<separateSpecs.length;i++){
			System.arraycopy(params, start, separateSpecs[i], 0, separateSpecs[i].length);
			start += separateSpecs[i].length;
		}
		for(int i=0;i<conditionalSpecs.length;i++){
			for(int j=0;j<conditionalSpecs[i].length;j++){
				System.arraycopy(params, start, conditionalSpecs[i][j], 0, conditionalSpecs[i][j].length);
				start += conditionalSpecs[i][j].length;
			}
		}

	}

	@Override
	public String getInstanceName() {
		return "LFSpecCond9C";
	}

	@Override
	public double getLogScoreFor(Sequence seq, int start) {
		
		ReferenceSequenceAnnotation data_anno =(ReferenceSequenceAnnotation)seq.getSequenceAnnotationByType("reference", 0);
		RVDSequence rvd_seq=(RVDSequence)data_anno.getReferenceSequence();
		
		int index_rj = rvd_seq.discreteVal(start-1);//current RVD
		
		int index_dj=rvd_seq.discreteValThirteen( index_rj );//13 current
		int index_zj=rvd_seq.discreteValTwelve( index_rj );//12 current
		
		int index_rjp = -1;//previous RVD
		int index_djp=-1;//13 previous
		int index_zjp=-1;//12 previous
		if(start-2 >= 0){
			index_rjp = rvd_seq.discreteVal(start-2);
			index_djp=rvd_seq.discreteValThirteen( index_rjp );
			index_zjp=rvd_seq.discreteValTwelve( index_rjp );
		}
		
		
		
		double score = specsThirteen[index_dj][seq.discreteVal(start)];
		
		if(separateMap[index_rj] >= 0){
			score += separateSpecs[separateMap[index_rj]][seq.discreteVal(start)];
		}
		
		if(start == 1 && firstPosMap[index_dj] >= 0){//first RVD and specific probability for that 13th AA
			
			score += firstPosSpecs[firstPosMap[index_dj]][seq.discreteVal(start)];
			
		}else if(index_zjp >= 0 && conditionalMapTwelve[index_zjp] >= 0 && conditionalMapRVD[index_rj] >= 0){
			
			score += conditionalSpecs[conditionalMapTwelve[index_zjp]][conditionalMapRVD[index_rj]][seq.discreteVal(start)];
			
		}
		
		return score;
		
	}

	@Override
	public boolean isInitialized() {
		return true;
	}

	@Override
	public String toString(NumberFormat nf) {
		
		StringBuffer sb = new StringBuffer();
		
		for(int i=0;i<specsThirteen.length;i++){
			sb.append(thirteen.getSymbol(0, i)+": "+Arrays.toString(specsThirteen[i])+"\n");
		}
		sb.append("\n");
		
		for(int i=0;i<separateMap.length;i++){
			if(separateMap[i] >= 0){
				sb.append(rvds.getSymbol(0, i)+": "+Arrays.toString(separateSpecs[separateMap[i]])+"\n");
			}
		}
		sb.append("\n");
		
		for(int i=0;i<firstPosMap.length;i++){
			if(firstPosMap[i] >= 0){
				sb.append(thirteen.getSymbol(0, i)+": "+Arrays.toString(firstPosSpecs[firstPosMap[i]])+"\n");
			}
		}
		sb.append("\n");
		
		for(int i=0;i<conditionalMapTwelve.length;i++){
			if(conditionalMapTwelve[i] >= 0){
				for(int j=0;j<conditionalMapRVD.length;j++){
					if(conditionalMapRVD[j] >= 0){
						sb.append(rvds.getSymbol(0, j)+"|"+twelve.getSymbol(0, i)+": "+Arrays.toString(conditionalSpecs[conditionalMapTwelve[i]][conditionalMapRVD[j]])+"\n");
					}
				}
			}
		}
	
		
		return sb.toString();
		
	}

	@Override
	public StringBuffer toXML() {
		StringBuffer xml = new StringBuffer();
		XMLParser.appendObjectWithTags(xml, conditionalMapRVD, "conditionalMapRVD");
		XMLParser.appendObjectWithTags(xml, conditionalMapTwelve, "conditionalMapTwelve");
		XMLParser.appendObjectWithTags(xml, conditionalSpecs, "conditionalSpecs");
		XMLParser.appendObjectWithTags(xml, firstPosMap, "firstPosMap");
		XMLParser.appendObjectWithTags(xml, firstPosSpecs, "firstPosSpecs");
		XMLParser.appendObjectWithTags(xml, offConditionalSpecs, "offConditionalSpecs");
		XMLParser.appendObjectWithTags(xml, offSeparateSpecs, "offSeparateSpecs");
		XMLParser.appendObjectWithTags(xml, offSpecsThirteen, "offSpecsThirteen");
		XMLParser.appendObjectWithTags(xml, separateMap, "separateMap");
		XMLParser.appendObjectWithTags(xml, separateSpecs, "separateSpecs");
		XMLParser.appendObjectWithTags(xml, specsThirteen, "specsThirteen");
		
		XMLParser.appendObjectWithTags(xml, twelve, "twelve");
		XMLParser.appendObjectWithTags(xml, thirteen, "thirteen");
		XMLParser.appendObjectWithTags(xml, rvds, "rvds");
		
		XMLParser.addTags(xml, "LFSpec9C");
		return xml;
	}

	@Override
	protected void fromXML(StringBuffer xml) throws NonParsableException {
		
		xml = XMLParser.extractForTag(xml, "LFSpec9C");
		conditionalMapRVD = (int[]) XMLParser.extractObjectForTags(xml, "conditionalMapRVD");
		conditionalMapTwelve = (int[]) XMLParser.extractObjectForTags(xml, "conditionalMapTwelve");
		conditionalSpecs = (double[][][]) XMLParser.extractObjectForTags(xml, "conditionalSpecs");
		firstPosMap = (int[]) XMLParser.extractObjectForTags(xml, "firstPosMap");
		firstPosSpecs = (double[][]) XMLParser.extractObjectForTags(xml, "firstPosSpecs");
		offConditionalSpecs = (int) XMLParser.extractObjectForTags(xml, "offConditionalSpecs");
		offSeparateSpecs = (int) XMLParser.extractObjectForTags(xml, "offSeparateSpecs");
		offSpecsThirteen = (int) XMLParser.extractObjectForTags(xml, "offSpecsThirteen");
		separateMap = (int[]) XMLParser.extractObjectForTags(xml, "separateMap");
		separateSpecs = (double[][]) XMLParser.extractObjectForTags(xml, "separateSpecs");
		specsThirteen = (double[][]) XMLParser.extractObjectForTags(xml, "specsThirteen");
		
		twelve = (AlphabetContainer) XMLParser.extractObjectForTags(xml, "twelve");
		thirteen = (AlphabetContainer) XMLParser.extractObjectForTags(xml, "thirteen");
		rvds = (AlphabetContainer) XMLParser.extractObjectForTags(xml, "rvds");
		
		alphabets = DNAAlphabetContainer.SINGLETON;
		length = 1;
	}

	public double[] getSpecs(RVDSequence rvd_seq, int start) {
		
		double[] specs = new double[(int) alphabets.getAlphabetLengthAt(0)];
		
		
		
		int index_rj = rvd_seq.discreteVal(start-1);//current RVD
		
		int index_dj=rvd_seq.discreteValThirteen( index_rj );//13 current
		int index_zj=rvd_seq.discreteValTwelve( index_rj );//12 current
		
		int index_rjp = -1;//previous RVD
		int index_djp=-1;//13 previous
		int index_zjp=-1;//12 previous
		if(start-2 >= 0){
			index_rjp = rvd_seq.discreteVal(start-2);
			index_djp=rvd_seq.discreteValThirteen( index_rjp );
			index_zjp=rvd_seq.discreteValTwelve( index_rjp );
		}
		
		
		for(int i=0;i<specs.length;i++){
		
		
			double score = specsThirteen[index_dj][i];
			
			if(separateMap[index_rj] >= 0){
				score += separateSpecs[separateMap[index_rj]][i];
			}
			
			if(start == 1 && firstPosMap[index_dj] >= 0){//first RVD and specific probability for that 13th AA
				
				score += firstPosSpecs[firstPosMap[index_dj]][i];
				
			}else if(index_zjp >= 0 && conditionalMapTwelve[index_zjp] >= 0 && conditionalMapRVD[index_rj] >= 0){
				
				score += conditionalSpecs[conditionalMapTwelve[index_zjp]][conditionalMapRVD[index_rj]][i];
				
			}
			
			specs[i] = score;
		
		}
		
		return specs;
		
		
		
	}

}
