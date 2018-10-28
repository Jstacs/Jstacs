package projects.xanthogenomes;

import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.HashSet;
import java.util.LinkedList;

import de.jstacs.data.DNADataSet;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.data.sequences.annotation.SimpleSequenceAnnotationParser;
import de.jstacs.sequenceScores.statisticalModels.differentiable.homogeneous.HomogeneousMMDiffSM;
import de.jstacs.sequenceScores.statisticalModels.trainable.hmm.AbstractHMM;
import de.jstacs.sequenceScores.statisticalModels.trainable.hmm.HMMFactory;
import de.jstacs.utils.ComparableElement;
import de.jstacs.utils.DoubleList;
import de.jstacs.utils.IntList;
import de.jstacs.utils.Pair;
import de.jstacs.utils.ToolBox;
import projects.xanthogenomes.tools.TALEPredictionTool;


public class NHMMerLoose {

	public static void main(String[] args) throws Exception {
		
		SimpleSequenceAnnotationParser parser = new SimpleSequenceAnnotationParser();
		
		DNADataSet ds = new DNADataSet( args[0],'>', parser );
		
		
		StringBuffer repeatConsensus = new StringBuffer();
		Pair<AbstractHMM, HomogeneousMMDiffSM> repeats = HMMFactory.parseProfileHMMFromHMMer( new InputStreamReader( TALEPredictionTool.class.getClassLoader().getResourceAsStream( "projects/xanthogenomes/data/repeats.hmm" ) ), repeatConsensus, null, null );
		
		StringBuffer startConsensus = new StringBuffer();
		LinkedList<Integer> startMatchStates = new LinkedList<Integer>();
		LinkedList<Integer> startSilentStates = new LinkedList<Integer>();
		Pair<AbstractHMM, HomogeneousMMDiffSM> start = HMMFactory.parseProfileHMMFromHMMer( new InputStreamReader( TALEPredictionTool.class.getClassLoader().getResourceAsStream("projects/xanthogenomes/data/starts.hmm")), startConsensus, startMatchStates, startSilentStates );
		
		StringBuffer endConsensus = new StringBuffer();
		LinkedList<Integer> endMatchStates = new LinkedList<Integer>();
		LinkedList<Integer> endSilentStates = new LinkedList<Integer>();
		Pair<AbstractHMM, HomogeneousMMDiffSM> end = HMMFactory.parseProfileHMMFromHMMer( new InputStreamReader( TALEPredictionTool.class.getClassLoader().getResourceAsStream("projects/xanthogenomes/data/ends.hmm")), endConsensus, endMatchStates, endSilentStates );
		
		
		PrintWriter wr = new PrintWriter(args[0]+"_stretch.fasta");
		for(int i=0;i<ds.getNumberOfElements();i++){
			Sequence seq = ds.getElementAt(i);
		//	System.out.println(i+":");
			int[][] res = run(repeats, start, end, repeatConsensus, startConsensus, endConsensus, startMatchStates, startSilentStates, endMatchStates, endSilentStates, seq);
			
			if(res.length > 0){

				int first = seq.getLength();
				int last =0;

				for(int j=0;j<res.length;j++){
					if(res[j][2] > 0){
						if(res[j][0] < first){
							first = res[j][0]; 
						}
						if(res[j][1] > last){
							last = res[j][1];
						}
					}else{
						int temp1 = seq.getLength()-res[j][1];
						int temp2 = seq.getLength()-res[j][0];
						if(temp1 < first){
							first = temp1;
						}
						if(temp2 > last){
							last = temp2;
						}
					}
				}


				if(first > 1000 && seq.getLength()-last > 1000){
					System.out.println(i+" "+first+" "+last+" "+seq.getLength());
					wr.println(parser.parseAnnotationToComment('>', seq.getAnnotation()));
					wr.println(seq.toString());
					wr.flush();
				}
			}
			
		}
		wr.close();
		
	}
	
	
	public static int[][] run(Pair<AbstractHMM,HomogeneousMMDiffSM> repeats, Pair<AbstractHMM,HomogeneousMMDiffSM>  start, Pair<AbstractHMM,HomogeneousMMDiffSM> end, 
			StringBuffer repeatConsensus, StringBuffer startConsensus, StringBuffer endConsensus, 
			LinkedList<Integer> startMatchStates, LinkedList<Integer> startSilentStates,
			LinkedList<Integer> endMatchStates, LinkedList<Integer> endSilentStates, Sequence seq) throws Exception {
		
				
		LinkedList<int[]> fwd = findRepeats( seq, repeats.getFirstElement(), repeats.getSecondElement(), repeatConsensus.toString() );
		LinkedList<int[]> rev = findRepeats( seq.reverseComplement(), repeats.getFirstElement(), repeats.getSecondElement(), repeatConsensus.toString() );
		
		int totalNum = fwd.size()+rev.size();
		
		double k=0;
		
		
		LinkedList<int[]> list = new LinkedList<int[]>();
		
		for(int i=0;i<fwd.size();i++,k++){
			int[] curr = fwd.get( i );
			//System.out.println("currfwd: "+Arrays.toString( curr ));
			int[] refinestartend = new int[]{curr[1],curr[2]};
		//	System.out.println("repeats: "+Arrays.toString(refinestartend));
			int[] nterm = getBestTerminus( seq, curr[1], curr[2], true, true, start.getFirstElement(), start.getSecondElement(), startConsensus.toString(), startMatchStates, startSilentStates );
			if(nterm != null){
				refinestartend[0] = nterm[0];
			}
			int[] cterm = getBestTerminus( seq, curr[1], curr[2], true, false, end.getFirstElement(), end.getSecondElement(), endConsensus.toString(), endMatchStates, endSilentStates );
			if(cterm != null){
				refinestartend[1] = cterm[1];
			}
			
			list.add(new int[]{refinestartend[0],refinestartend[1],+1});
			
		}
		
		for(int i=0;i<rev.size();i++, k++){
			int[] curr = rev.get( i );
			//System.out.println("currrev: "+Arrays.toString( curr ));
			int[] refinestartend = new int[]{curr[1],curr[2]};
			//System.out.println("repeats: "+Arrays.toString(refinestartend));
			int[] nterm = getBestTerminus( seq, curr[1], curr[2], false, true, start.getFirstElement(), start.getSecondElement(), startConsensus.toString(), startMatchStates, startSilentStates );
			if(nterm != null){
				refinestartend[0] = nterm[0];
			}
			int[] cterm = getBestTerminus( seq, curr[1], curr[2], false, false, end.getFirstElement(), end.getSecondElement(), endConsensus.toString(), endMatchStates, endSilentStates );
			if(cterm != null){
				refinestartend[1] = cterm[1];
			}
			
			list.add(new int[]{refinestartend[0],refinestartend[1],-1});
		}
		
		LinkedList<int[]> toRemove = new LinkedList<int[]>();
		
		for(int i=1;i<list.size();i++){
			int[] temp = list.get(i-1);
			int[] temp2 = list.get(i);
			if(temp[2] == temp2[2]){
				if(temp[0] >= temp2[0] && temp[1] <= temp2[1]){
					toRemove.add(temp);
				}else if(temp[0] <= temp2[0] && temp[1] >= temp2[1]){
					toRemove.add(temp2);
				}
			}
		}
		list.removeAll(toRemove);
		
		return list.toArray( new int[0][] );
	}
	

	public static int[] getBestTerminus(Sequence seq, int start, int end, boolean fwd, boolean isStart, AbstractHMM hmm, HomogeneousMMDiffSM hom, String consensus, LinkedList<Integer> matchStates, LinkedList<Integer> silentStates) throws Exception {
		
		int numLay = consensus.length();//TODO
		numLay = (int)Math.round(numLay*1.1);
		
		int w = numLay;
		
		double t = consensus.length()*Math.log( 1.1 );
		
		
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
	
	public static LinkedList<int[]> findRepeats(Sequence seq, AbstractHMM hmm, HomogeneousMMDiffSM hom, String consensus) throws Exception {
		
		int totalLength = seq.getLength();
		
		int numLay = consensus.length();
		numLay = (int)Math.round(numLay*1.1);
		
		int w = numLay;
		
		int frag = 10;
		HashSet<String> parts = new HashSet<String>();
		for(int i=0;i<consensus.length()/frag;i++){
			parts.add( consensus.substring( i*frag, (i+1 )*frag).toUpperCase() );
		}
		
		
		double t = consensus.length()*Math.log( 1.1 );
		
		
		
		LinkedList<int[]> found = new LinkedList<int[]>();
		
		double l=0;
		
		
			
			double[] vals = new double[seq.getLength()-w+1];
			
			int num = -1;
			
			for(int j=0;j<seq.getLength()-w+1;j++,l++){
				
				
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
				
				
				if(num > parts.size()/4){
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
				//System.out.println("max:"+max+" > "+t);
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
				if(els.length > 1){
					for(int j=1;j<els.length;j++){
						//System.out.println(els[j].getWeight());
						if(els[j].getWeight()-1000 > els[j-1].getWeight() || j == els.length-1){
							end = els[j-1].getWeight()+consensus.length();
							//System.out.println(start+" "+end );
							found.add( new int[]{0,start,end} );
							start = els[j].getWeight();
							end = els[j].getWeight()+consensus.length();
						}
					}
				}else{
					found.add(new int[]{0,start,end});
				}
			}
		
		//System.out.println("time: "+time.getElapsedTime());
		return found;
	}
	
}
