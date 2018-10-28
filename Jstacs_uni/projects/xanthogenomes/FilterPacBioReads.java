package projects.xanthogenomes;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.io.Reader;
import java.util.HashSet;

import de.jstacs.data.alphabets.DNAAlphabetContainer;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.sequenceScores.statisticalModels.differentiable.homogeneous.HomogeneousMMDiffSM;
import de.jstacs.sequenceScores.statisticalModels.trainable.hmm.AbstractHMM;
import de.jstacs.sequenceScores.statisticalModels.trainable.hmm.HMMFactory;
import de.jstacs.utils.Pair;
import projects.xanthogenomes.tools.TALEPredictionTool;

public class FilterPacBioReads {

	public static void main(String[] args) throws Exception {
		
		
		Reader repeatHMMer =  new InputStreamReader( TALEPredictionTool.class.getClassLoader().getResourceAsStream( "projects/xanthogenomes/data/repeats.hmm" ) );
		
		StringBuffer consensus = new StringBuffer();
		Pair<AbstractHMM, HomogeneousMMDiffSM> repeats = HMMFactory.parseProfileHMMFromHMMer( repeatHMMer, consensus, null, null );
		
		int frag = 10;
		HashSet<String> parts = new HashSet<String>();
		for(int i=0;i<consensus.length()/frag;i++){
			parts.add( consensus.substring( i*frag, (i+1 )*frag).toUpperCase() );
		}
		
		BufferedReader read = new BufferedReader(new FileReader(args[0]));
		
		PrintWriter no = new PrintWriter(args[0]+"_norepeats.fastq");
		PrintWriter re = new PrintWriter(args[0]+"_repeats.fastq");
		
		String str = null;
		int i=0;
		
		String head = null;
		String seq = null;
		String qual = null;
		boolean use = false;
		
		
		
		while( (str = read.readLine()) != null ){
			if(i%4==0){
				head = str;
			}else if( (i-1)%4==0 ){
				seq = str;
				use = ! findRepeats(Sequence.create(DNAAlphabetContainer.SINGLETON, seq), repeats.getFirstElement(), repeats.getSecondElement(), parts, frag, consensus.length());
				use &= ! findRepeats(Sequence.create(DNAAlphabetContainer.SINGLETON, seq).reverseComplement(), repeats.getFirstElement(), repeats.getSecondElement(), parts, frag, consensus.length());
			}else if( (i-3)%4==0 ){
				qual = str;
				PrintWriter temp = null;
				if(use){
					temp = no;
				}else{
					temp = re;
				}
				temp.println(head);
				temp.println(seq);
				temp.println("+");
				temp.println(qual);
			}
			i++;
		}
		
		no.close();
		re.close();
		read.close();

	}




	public static boolean findRepeats(Sequence seq, AbstractHMM hmm, HomogeneousMMDiffSM hom, HashSet<String> parts, int frag, int consensusLength) throws Exception {

		int numLay = consensusLength;
		numLay = (int)Math.round(numLay*1.1);

		int w = numLay;




		double t = consensusLength*Math.log( 1.2 );

		double maxVal=Double.NEGATIVE_INFINITY;

		int num = -1;

		for(int j=0;j<seq.getLength()-w+1;j++){


			Sequence sub = seq.getSubSequence( j, w );


			if(num == -1){
				num = 0;
				String substr = sub.toString();
				String[] parts2 = parts.toArray( new String[0] );
				for(int k=0;k<parts2.length;k++){
					if(substr.indexOf( parts2[k] ) > -1){
						num++;
					}
				}
			}


			if(num > parts.size()/3){
				double fg = hmm.getLogProbFor( sub );
				double bg = hom.getLogProbFor( sub );
				double rat = fg-bg;
				if(rat > maxVal){
					maxVal = rat;
				}
			}


			if(j<seq.getLength()-w){
				String substr1 = seq.toString( j, j+frag );
				String substr2 = seq.toString( j+w-frag+1,j+w+1 );
				if(parts.contains( substr1 )){
					num--;
				}
				if(parts.contains( substr2 )){
					num++;
				}
			}
		}
		//System.out.println(maxVal+" > "+t);
		if(maxVal > t){
			return true;
		}else{
			return false;
		}
	}

}
