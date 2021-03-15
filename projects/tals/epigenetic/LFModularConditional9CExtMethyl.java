package projects.tals.epigenetic;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.StringReader;
import java.util.Arrays;

import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.WrongAlphabetException;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.data.sequences.annotation.ReferenceSequenceAnnotation;
import de.jstacs.data.sequences.annotation.SequenceAnnotation;
import de.jstacs.io.FileManager;
import projects.tals.RVDSequence;
import projects.tals.linear.LFModularConditional9C;

public class LFModularConditional9CExtMethyl extends LFModularConditional9C {
	private double[] methylSpecsThirteen; 
	private int[] separateMap;
	private double[] methylSeparateSpecs;
	
	public LFModularConditional9CExtMethyl(StringBuffer xml,AlphabetContainer thirteen,AlphabetContainer rvds,String[] separateRVDs) throws Exception {
		super(xml);
		this.methylSpecsThirteen=new double[(int)thirteen.getAlphabetLengthAt(0)];
		this.methylSeparateSpecs = new double[separateRVDs.length];
		StringBuffer SF_13Specs_MethylC=FileManager.readInputStream( QuickTBSPredictionToolMethylAccessibilityAnnotation_fai.class.getClassLoader().getResourceAsStream( "projects/tals/epigenetic/PrediTALE13Specs_MethylC.csv"));
		BufferedReader br = new BufferedReader(new BufferedReader(new StringReader(SF_13Specs_MethylC.toString())));
		String line="";
		String[] splitLine;
		while ((line = br.readLine()) != null){
			splitLine=line.split(",");
			methylSpecsThirteen[(int) thirteen.getCode(0, splitLine[0])]=Double.parseDouble(splitLine[1]);
		}
		br.close();
		StringBuffer SF_comRVDSpecs_MethylC=FileManager.readInputStream( QuickTBSPredictionToolMethylAccessibilityAnnotation_fai.class.getClassLoader().getResourceAsStream( "projects/tals/epigenetic/PrediTALEcomRVDSpecs_MethylC.csv"));
		BufferedReader br2 = new BufferedReader(new BufferedReader(new StringReader(SF_comRVDSpecs_MethylC.toString())));
		
		separateMap = new int[(int) rvds.getAlphabetLengthAt(0)];
		Arrays.fill(separateMap, -1);
		int k=0;
		for(String sepRVD : separateRVDs){
			separateMap[ (int) rvds.getCode(0, sepRVD) ] = k;
			k++;
		}

		line="";
		while ((line = br2.readLine()) != null){
			
			splitLine=line.split(",");
			
			if(Arrays.asList(separateRVDs).contains(splitLine[0])){
				methylSeparateSpecs[separateMap[(int) rvds.getCode(0, splitLine[0])]]=Double.parseDouble(splitLine[1]);
			}
			
		}
		br2.close();
	}
	
	public double[][] toPWM(RVDSequence rvds){
		double[][] pwm = new double[rvds.getLength()+1][];
		double[] methylC_specs =new double[rvds.getLength()+1];
		pwm[0] = lf0.getSpecs(rvds);
		
		for(int i=1;i<pwm.length;i++){
			double[] specs=specificity.getSpecs(rvds,i);
			int index_rj = rvds.discreteVal(i-1);//current RVD
			int index_dj=rvds.discreteValThirteen( index_rj );//13 current
			methylC_specs[i] = methylSpecsThirteen[index_dj];
			if(separateMap[index_rj] >= 0){
				methylC_specs[i] += methylSeparateSpecs[separateMap[index_rj]];
			}
			pwm[i]=new double[specs.length+1];
			System.arraycopy(specs, 0, pwm[i], 0, specs.length);
			pwm[i][specs.length]=methylC_specs[i];
			
			double pos = (position == null ? 1.0 : position.getLogScoreFor(rvds.getLength()+1, i));
			
			for(int j=0;j<pwm[i].length;j++){
				pwm[i][j] *= pos;
			}
		}
		
		for(int i=0;i<pwm.length;i++){
			for(int j=0;j<pwm[i].length;j++){
				pwm[i][j] /= pwm.length;
			}
		}
		//System.out.println(rvds.toString());
		//for (int i = 0; i < pwm.length; i++) {
			//System.out.println(Arrays.toString(pwm[i]));
		//}
		return pwm;
		
	}
	
	public double getLogScoreFor(Sequence seq, int start) {
		MethylationSequenceAnnotation methylAnno = (MethylationSequenceAnnotation)seq.getSequenceAnnotationByType("methylationprofil", 0);
		Methylationprofil MP=methylAnno.getMethylationprofile();
	
		ReferenceSequenceAnnotation data_anno =(ReferenceSequenceAnnotation)seq.getSequenceAnnotationByType("reference", 0);
		RVDSequence rvd_seq=(RVDSequence)data_anno.getReferenceSequence();
		
		String mask = null;
		SequenceAnnotation mann = seq.getSequenceAnnotationByType("mask", 0);
		if(mann != null){
			mask = mann.getIdentifier();
		}
		SequenceAnnotation ann = seq.getSequenceAnnotationByType("intgroup", 0);
		int group = 0;
		if(ann != null){
			String gs = ann.getIdentifier();
			group = Integer.parseInt( gs );
		}
		
		double score = 0.0;
		if(mask == null || mask.charAt(start) == 'O'){
			score += lf0.getLogScoreFor(seq, start);
		}
		
		for(int i=start+1;i<seq.getLength();i++){
			if(mask == null || mask.charAt(i) == 'O'){
				
				int index_rj = rvd_seq.discreteVal(i-1);//current RVD
				int index_dj=rvd_seq.discreteValThirteen( index_rj );//13 current
				double spec = (1-MP.getMethylPropAtPos(i))*specificity.getLogScoreFor(seq, i);

				try {
					if((MP.getMethylPropAtPos(i)>0)&(seq.discreteVal( i )==seq.getAlphabetContainer().getCode(0, "C"))){
						
						double methylSpecificity=methylSpecsThirteen[index_dj];
						if(separateMap[index_rj] >= 0){
							methylSpecificity += methylSeparateSpecs[separateMap[index_rj]];
						}
						spec+=(MP.getMethylPropAtPos(i))*methylSpecificity;
					}
				} catch (WrongAlphabetException e) {
					e.printStackTrace();
				}
				
				double pos = (position == null ? 1.0 : position.getLogScoreFor(seq, i));
				
				score += spec*pos;
				
			}
		}
		
		score *= Math.exp( a[group] );
		score += b[group];
		
		return score;
	}
}
