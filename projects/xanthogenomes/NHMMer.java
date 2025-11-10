package projects.xanthogenomes;

import java.io.FileReader;
import java.io.Reader;
import java.util.Arrays;
import java.util.HashSet;
import java.util.LinkedList;

import de.jstacs.data.DNADataSet;
import de.jstacs.data.DataSet;
import de.jstacs.data.WrongAlphabetException;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.data.sequences.WrongSequenceTypeException;
import de.jstacs.sequenceScores.statisticalModels.differentiable.homogeneous.HomogeneousMMDiffSM;
import de.jstacs.sequenceScores.statisticalModels.trainable.hmm.AbstractHMM;
import de.jstacs.sequenceScores.statisticalModels.trainable.hmm.HMMFactory;
import de.jstacs.tools.ProgressUpdater;
import de.jstacs.utils.ComparableElement;
import de.jstacs.utils.DoubleList;
import de.jstacs.utils.IntList;
import de.jstacs.utils.Pair;
import de.jstacs.utils.ToolBox;
import projects.xanthogenomes.Tools.Translator;


public class NHMMer {

	public static void main(String[] args) throws Exception {
		DNADataSet ds = new DNADataSet( args[3] );
		
		int[][] res = run(new FileReader(args[0]), new FileReader( args[1]), new FileReader( args[2] ), ds, new ProgressUpdater(), false);
		
		for(int i=0;i<res.length;i++){
			System.out.println(Arrays.toString( res[i] ));
		}
		
		
	}
	
	
	public static int[][] run(Reader repeatHMMer, Reader startHMMer, Reader endHMMer, DataSet ds, ProgressUpdater progress, boolean sensitive) throws Exception {
		StringBuffer repeatConsensus = new StringBuffer();
		Pair<AbstractHMM, HomogeneousMMDiffSM> repeats = HMMFactory.parseProfileHMMFromHMMer( repeatHMMer, repeatConsensus, null, null );
		
		StringBuffer startConsensus = new StringBuffer();
		LinkedList<Integer> startMatchStates = new LinkedList<Integer>();
		LinkedList<Integer> startSilentStates = new LinkedList<Integer>();
		Pair<AbstractHMM, HomogeneousMMDiffSM> start = HMMFactory.parseProfileHMMFromHMMer( startHMMer, startConsensus, startMatchStates, startSilentStates );
		
		StringBuffer endConsensus = new StringBuffer();
		LinkedList<Integer> endMatchStates = new LinkedList<Integer>();
		LinkedList<Integer> endSilentStates = new LinkedList<Integer>();
		Pair<AbstractHMM, HomogeneousMMDiffSM> end = HMMFactory.parseProfileHMMFromHMMer( endHMMer, endConsensus, endMatchStates, endSilentStates );
		
		progress.setLast( 4.0 );
		
		LinkedList<int[]> fwd = findRepeats( ds, repeats.getFirstElement(), repeats.getSecondElement(), repeatConsensus.toString(), progress, 0.0, sensitive );
		LinkedList<int[]> rev = findRepeats( ds.getReverseComplementaryDataSet(), repeats.getFirstElement(), repeats.getSecondElement(), repeatConsensus.toString(), progress, 1.0, sensitive );
		
		int totalNum = fwd.size()+rev.size();
		
		double k=0;
		
		
		LinkedList<int[]> list = new LinkedList<int[]>();
		
		for(int i=0;i<fwd.size();i++,k++){
			int[] curr = fwd.get( i );
			//System.out.println("curr: "+Arrays.toString( curr ));
			int[] refinestartend = new int[]{curr[1],curr[2]};
		//	System.out.println("repeats: "+Arrays.toString(refinestartend));
			int[] nterm = getBestTerminus( ds, curr[0], curr[1], curr[2], true, true, start.getFirstElement(), start.getSecondElement(), startConsensus.toString(), startMatchStates, startSilentStates );
			if(nterm != null){
				refinestartend[0] = nterm[0];
			}
			int[] cterm = getBestTerminus( ds, curr[0], curr[1], curr[2], true, false, end.getFirstElement(), end.getSecondElement(), endConsensus.toString(), endMatchStates, endSilentStates );
			if(cterm != null){
				refinestartend[1] = cterm[1];
			}
			//System.out.println("terms: "+Arrays.toString( nterm )+" "+Arrays.toString( cterm ));
			int[] reg = refine(refinestartend[0], refinestartend[1],ds.getElementAt( curr[0] ));
			//System.out.println("refined: "+Arrays.toString( reg ));
			
			int mRNALength = Math.max( refinestartend[1], reg[1] ) - Math.min( refinestartend[0], reg[0] );
			int cdsLength = reg[1] - reg[0];
			
			if(mRNALength - startConsensus.length()/3 > cdsLength || mRNALength - endConsensus.length()/3 > cdsLength){
				list.add( new int[]{curr[0], reg[0],reg[1],+1, Math.min( refinestartend[0], reg[0] ), Math.max( refinestartend[1], reg[1] ), 1} );
			}else{
				list.add( new int[]{curr[0], reg[0],reg[1],+1, reg[0], reg[1], 0} );
			}
			
			progress.setCurrent( 2.0*k/totalNum + 2.0 );
			
		}
		
		for(int i=0;i<rev.size();i++, k++){
			int[] curr = rev.get( i );
			//System.out.println("curr: "+Arrays.toString( curr ));
			int[] refinestartend = new int[]{curr[1],curr[2]};
			//System.out.println("repeats: "+Arrays.toString(refinestartend));
			int[] nterm = getBestTerminus( ds, curr[0], curr[1], curr[2], false, true, start.getFirstElement(), start.getSecondElement(), startConsensus.toString(), startMatchStates, startSilentStates );
			if(nterm != null){
				refinestartend[0] = nterm[0];
			}
			int[] cterm = getBestTerminus( ds, curr[0], curr[1], curr[2], false, false, end.getFirstElement(), end.getSecondElement(), endConsensus.toString(), endMatchStates, endSilentStates );
			if(cterm != null){
				refinestartend[1] = cterm[1];
			}
			//System.out.println("terms: "+Arrays.toString( nterm )+" "+Arrays.toString( cterm ));
			int[] reg = refine( refinestartend[0], refinestartend[1], ds.getElementAt( curr[0] ).reverseComplement() );
			//System.out.println("refined: "+Arrays.toString( reg ));
			int el = ds.getElementAt( curr[0] ).getLength();
			
			int mRNALength = Math.max( refinestartend[1], reg[1] ) - Math.min( refinestartend[0], reg[0] );
			int cdsLength = reg[1] - reg[0];
			
			if(mRNALength - startConsensus.length()/3 > cdsLength || mRNALength - endConsensus.length()/3 > cdsLength){
				list.add( new int[]{curr[0],el-reg[1],el-reg[0],-1,el-Math.max( refinestartend[1], reg[1] ),el-Math.min( refinestartend[0], reg[0] ),1} );
			}else{
				list.add( new int[]{curr[0],el-reg[1],el-reg[0],-1,el-reg[1],el-reg[0],0} );
			}
			progress.setCurrent( 2.0*k/totalNum + 2.0 );
		}
		
		LinkedList<int[]> toRemove = new LinkedList<int[]>();
		
		for(int i=1;i<list.size();i++){
			int[] temp = list.get(i-1);
			int[] temp2 = list.get(i);
			if(temp[0] == temp2[0] && temp[3] == temp2[3]){
				if(temp[1] >= temp2[1] && temp[2] <= temp2[2]){
					toRemove.add(temp);
				}else if(temp[1] <= temp2[1] && temp[2] >= temp2[2]){
					toRemove.add(temp2);
				}
			}
		}
		list.removeAll(toRemove);
		
		return list.toArray( new int[0][] );
	}
	
	private static int[] refine( int start, int end, Sequence original ) throws WrongAlphabetException, WrongSequenceTypeException {
		Sequence sequence = original.getSubSequence( start, end-start+1 );
		int maxframe = -1;
		int maxorf = 0;
		int maxlen = 0;
		for(int i=0;i<3;i++){
			Sequence prot = Translator.DEFAULT.translate( sequence, i );
			//System.out.println("frame "+i+" "+prot.toString());
			String[] parts = prot.toString().split( "\\*" );
			for(int j=0;j<parts.length;j++){
				if(parts[j].length() > maxlen){
					maxlen = parts[j].length();
					maxframe = i;
					maxorf = j;
				}
			}
		}
		//System.out.println("frame: "+maxframe+" "+maxorf+" "+maxlen);
		Sequence prot = Translator.DEFAULT.translate( sequence, maxframe );
		//System.out.println("prot: "+prot);
		String[] parts = prot.toString().split( "\\*" );
		int off = maxframe;
		for(int i=0;i<maxorf;i++){
			off += (parts[i].length()+1)*3;
		}
		int len = (parts[maxorf].length()+1)*3;
	//	System.out.println("before "+off+" "+len);
		if(maxorf > 0){
			Sequence sub = original.getSubSequence( start+off, len );
			for(int i=0;i<sub.getLength()-2;i+=3){
				String codon = sub.toString( i, i+3 );
				if("ATG".equals( codon )){
					off += i;
					len -= i;
					break;
				}
			}
		}else{
			while( start+off-3 >= 0 && !"ATG".equals( original.toString( start+off, start+off+3 ) ) ){
				off -= 3;
				len += 3;
			}
		}
		//System.out.println("after: "+off+" "+len);
		while( start+off+len+3<original.getLength() && !"*".equals( Translator.DEFAULT.translate( original.getSubSequence( start+off+len-3, 3 ), 0 ).toString() ) ){
			len += 3;
		}
		
		while( start+off+len+3>original.getLength() ) {
			len -= 3;
		}
		
		//Sequence sub = original.getSubSequence( start+off, len );
		//System.out.println(sub.getLength()+" ");
		//System.out.println(Translator.DEFAULT.translate( sub, 0 ));
		
		return new int[]{start+off,start+off+len};
	}

	public static int[] getBestTerminus(DataSet ds, int id, int start, int end, boolean fwd, boolean isStart, AbstractHMM hmm, HomogeneousMMDiffSM hom, String consensus, LinkedList<Integer> matchStates, LinkedList<Integer> silentStates) throws Exception {
		
		int numLay = consensus.length();//TODO
		numLay = (int)Math.round(numLay*1.1);
		
		int w = numLay;
		
		double t = consensus.length()*Math.log( 1.5 );
		
		Sequence seq = ds.getElementAt( id );
		
		if(!fwd){
			seq = seq.reverseComplement();
		}
		
		DoubleList scores = new DoubleList();
		IntList positions = new IntList();
		
		if(isStart){
			int numLower = 0;
			for(int i=start-w+(int)Math.round( 0.1*consensus.length() );i>=Math.max( 0, start-w-200);i-=5){
				//Sequence sub = seq.getSubSequence( i, w );
				double fg = hmm.getLogProbFor( seq, i, i+w-1 );
				double bg = hom.getLogProbFor( seq, i, i+w-1 );
				double rat = fg-bg;
			//	System.out.println(i+" "+rat);
				if(scores.length() > 0 && rat < scores.get( scores.length()-1 )){
					numLower++;
				}else{
					numLower = 0;
				}
				if(numLower > 10){
					break;
				}
				scores.add( rat );
				positions.add( i );
			}
		}else{
			int numLower = 0;
			for(int i=end-(int)Math.round( 0.1*consensus.length() );i<Math.min( end+200, seq.getLength()-w+1 );i+=5){
				//Sequence sub = seq.getSubSequence( i, w );
				double fg = hmm.getLogProbFor( seq,i, i+w-1 );
				double bg = hom.getLogProbFor( seq, i, i+w-1 );
				double rat = fg-bg;
				if(scores.length() > 0 && rat < scores.get( scores.length()-1 )){
					numLower++;
				}else{
					numLower = 0;
				}
				if(numLower > 10){
					break;
				}
				//System.out.println(i+" "+rat+" "+numLower);
				scores.add( rat );
				positions.add( i );
			}
		}

		if(positions.length()==0){
			return null;
		}
		
		int idx = ToolBox.getMaxIndex( scores.toArray() );
		//System.out.println("poslen: "+positions.length()+", idx: "+idx+" "+isStart);
		//System.out.println("end: "+end+" conslen: "+consensus.length()+" "+seq.getLength()+" "+w);
		
		/*if(scores.get( idx )<t){
			return null;
		}else{*/


			int[] region = new int[2];

			Pair<IntList,Double> vit = hmm.getViterbiPathFor( positions.get( idx ), positions.get( idx )+w-1, seq );
			IntList states = vit.getFirstElement();
			double[] count = new double[states.length()];
			for(int i=0;i<states.length();i++){
				if(matchStates.contains( states.get( i ) ) ){
					count[i] = i > 0 ? count[i-1]+1 : 1;
				}else{
					count[i] = i > 0 && count[i-1]>0 ? count[i-1]-1 : 0;
				}
			}
			int endIdx = ToolBox.getMaxIndex( count );
			int startIdx = endIdx;
			for(;startIdx>=0 && count[startIdx]>0;startIdx--);
			//System.out.println(states);
			//System.out.println(Arrays.toString( count ));
			//System.out.println(startIdx+" "+endIdx);
			
			int offStart = 0;
			for(int i=0;i<startIdx;i++){
				if(!silentStates.contains( states.get( i )) ){
					offStart++;
				}
			}
			int offEnd = 0;
			for(int i=states.length()-1;i>endIdx;i--){
				if(!silentStates.contains( states.get( i )) ){
					offEnd++;
				}
			}
			//System.out.println("offs: "+offStart+" "+offEnd);
			region[0] = positions.get( idx )+offStart;
			region[1] = positions.get( idx )+w-offEnd;

			return region;
		//}
		
	}
	
	public static LinkedList<int[]> findRepeats(DataSet ds, AbstractHMM hmm, HomogeneousMMDiffSM hom, String consensus, ProgressUpdater progress, double progressOffset, boolean sensitive) throws Exception {
		
		int totalLength = 0;
		for(int i=0;i<ds.getNumberOfElements();i++){
			totalLength += ds.getElementAt( i ).getLength();
		}
		
		int numLay = consensus.length();
		numLay = (int)Math.round(numLay*1.1);
		
		int w = numLay;
		
		int frag = 10;
		if(sensitive) {
			frag = 5;
		}
		HashSet<String> parts = new HashSet<String>();
		for(int i=0;i<consensus.length()/frag;i++){
			parts.add( consensus.substring( i*frag, (i+1 )*frag).toUpperCase() );
		}
		
		
		double t = consensus.length()*Math.log( 1.3 );
		
		
		
		LinkedList<int[]> found = new LinkedList<int[]>();
		
		double l=0;
		
		for(int i=0;i<ds.getNumberOfElements();i++){
			Sequence seq = ds.getElementAt( i );

			if(seq.getLength()>=w) {

				double[] vals = new double[seq.getLength()-w+1];

				int num = -1;

				for(int j=0;j<seq.getLength()-w+1;j++,l++){

					if(j%1000 == 0){
						progress.setCurrent( l/totalLength + progressOffset );
					}

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

					if(num > parts.size()/2){
						double fg = hmm.getLogProbFor( sub );
						double bg = hom.getLogProbFor( sub );
						double rat = fg-bg;
						vals[j] = rat;
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
				//System.out.println("scan finished");
				LinkedList<ComparableElement<Double,Integer>> list = new LinkedList<ComparableElement<Double,Integer>>();
				while(true){
					int maxIdx = ToolBox.getMaxIndex( vals );
					double max = vals[maxIdx];
					for(int j=Math.max( 0, maxIdx-w/2);j<maxIdx+w/2 && j< vals.length; j++){
						vals[j] = 0;
					}
					if(max > t){
						list.add( new ComparableElement<Double,Integer>(max,maxIdx) );
					}else{
						break;
					}

				}

				ComparableElement<Double, Integer>[] els = list.toArray( new ComparableElement[0] );
				Arrays.sort( els );
				if(els.length > 0){
					int start = els[0].getWeight();
					int end = els[0].getWeight()+consensus.length();
					for(int j=1;j<els.length;j++){
						//System.out.println(els[j].getWeight());
						if(els[j].getWeight()-500 > els[j-1].getWeight() || j == els.length-1){
							end = els[j-1].getWeight()+consensus.length();
							//System.out.println(start+" "+end );
							found.add( new int[]{i,start,end} );
							start = els[j].getWeight();
							end = els[j].getWeight()+consensus.length();
						}
					}
				}
			}
		}
		
		//System.out.println("time: "+time.getElapsedTime());
		return found;
	}
	
}
