package projects.quickscan;

import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.Random;

import javax.naming.OperationNotSupportedException;

import de.jstacs.classifiers.differentiableSequenceScoreBased.gendismix.GenDisMixClassifier;
import de.jstacs.data.WrongAlphabetException;
import de.jstacs.data.alphabets.DiscreteAlphabet;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.data.sequences.WrongSequenceTypeException;
import de.jstacs.io.FileManager;
import de.jstacs.sequenceScores.QuickScanningSequenceScore;
import de.jstacs.utils.IntList;
import de.jstacs.utils.LargeSequenceReader;
import de.jstacs.utils.Pair;
import de.jstacs.utils.ToolBox;
import projects.dimont.ThresholdedStrandChIPper;

public class QuickMotifProfileTool {

	
	private static class ScoreRegion{
		
		private float[][] scores;
		private int index;
		private String id;
		private int offset;
		
		public ScoreRegion(float[][] scores, int index, String id, int offset) {
			super();
			this.scores = scores;
			this.index = index;
			this.id = id;
			this.offset = offset;
		}
		
		
		
	}
	
	private static class Lock{
		private int numFinished;
		private boolean[] isFinished;
		
		public Lock(int num){
			numFinished=num;
			isFinished = new boolean[num];
		}
		
		public synchronized void decrement(){
			numFinished--;
		}
		
		public synchronized boolean finished(){
			return numFinished==0;
		}
		
		public synchronized void setFinished(int num, boolean val){
			if(isFinished[num] == val){
				throw new RuntimeException();
			}else{
				if(isFinished[num] && !val){
					numFinished--;
				}else{
					numFinished++;
				}
				isFinished[num] = val;
			}
		}
		
		public synchronized int getFinished(){
			for(int i=0;i<isFinished.length;i++){
				if(isFinished[i]){
					return i;
				}
			}
			return -1;
		}
		
		public boolean allFinished(){
			boolean all = true;
			for(int i=0;all&&i<isFinished.length;i++){
				all &= isFinished[i];
			}
			return all;
		}
	}
	
	
	
	private int[] pow1;
	private int[] pow2;
	
	
	public static void main(String[] args) throws Exception{
		new QuickMotifProfileTool().run(args[0], args[1], Integer.parseInt(args[2]));
	}
	
	public void run(String slimfile, String genome, int numThreads) throws Exception {
		
		GenDisMixClassifier cl = new GenDisMixClassifier(FileManager.readFile(slimfile));
		ThresholdedStrandChIPper fg =(ThresholdedStrandChIPper) cl.getDifferentiableSequenceScore(0);
		QuickScanningSequenceScore lslim = (QuickScanningSequenceScore) fg.getFunction(0);
		
		BufferedReader read = new BufferedReader(new FileReader(genome));
		StringBuffer lastHeader = new StringBuffer();
		

		int a = (int)lslim.getAlphabetContainer().getAlphabetLengthAt(0);
		int prefK = 8;
		int restK = lslim.getLength()-prefK;
		int idxK = 4;
		
		pow1 = new int[prefK];
		pow1[pow1.length-1]=1;
		for( int i = pow1.length-2; i >= 0; i-- ) {
			pow1[i] = pow1[i+1]*a;
		}
		
		pow2 = new int[restK];
		pow2[pow2.length-1]=1;
		for( int i = pow2.length-2; i >= 0; i-- ) {
			pow2[i] = pow2[i+1]*a;
		}
		int maxValue = pow2[0]*a;
		
		long time = System.currentTimeMillis();
		
		IntList[] lists = new IntList[(int) Math.pow(a, prefK)];
		for(int i=0;i<lists.length;i++){
			lists[i] = new IntList(256);
		}
		
		Pair<IntList,ArrayList<Sequence>> pair = null;

		while( (pair = LargeSequenceReader.readNextSequences(read, lastHeader, lslim.getLength())) != null ){
			ArrayList<Sequence> seqs = pair.getSecondElement();
			Iterator<Sequence> it = seqs.iterator();
			int itIdx = 0;
			
			
			while( it.hasNext() ) {
				Sequence seq = it.next();
				
				String id = seq.getSequenceAnnotationByType("id", 0).getIdentifier().trim();
				System.err.println(id);
				
				itIdx++;
				//System.out.println(id);
				for(int d=0;d<2;d++){
					
					int idx1 = 0;
					for(int i=0;i<pow1.length-1;i++){
						idx1 += pow1[i+1]*seq.discreteVal(i);
					}
					int idx2 = 0;
					for(int i=0;i<pow2.length-1;i++){
						idx2 += pow2[i+1]*seq.discreteVal(i+prefK);
					}
					
					for(int j=0;j<seq.getLength()-lslim.getLength()+1;j++){
						idx1 = (idx1%pow1[0])*4 + seq.discreteVal(j+prefK-1);
						idx2 = (idx2%pow2[0])*4 + seq.discreteVal(j+prefK+restK-1);
						
						lists[idx1].add(idx2);
					}
					seq = seq.reverseComplement();
				}
				
			}
		}
		
		read.close();
		
		
		Random r = new Random();
		IntList[] assign = new IntList[numThreads];
		for(int i=0;i<assign.length;i++){
			assign[i] = new IntList();
		}
		for(int i=0;i<lists.length;i++){
			assign[r.nextInt(numThreads)].add(i);
		}
		
		
		
		
		int[][] acgt = new int[lists.length][(int) Math.pow(a, idxK)];
		
		int numPer2 = lists.length/numThreads;
		Lock lock3 = new Lock(numThreads);
		
		for(int i=0;i<numThreads;i++){
			final int l = i;
			new Thread(() -> {
				
				System.err.println("thread: "+l);
				for(int m=0;m<assign[l].length();m++){
					int k = assign[l].get(m);
					lists[k].sort();
					lists[k] = condense(lists[k]);
					
					for(int j=0;j<lists[k].length();j++){
						acgt[k][ lists[k].get(j)/pow2[idxK-1] ]++;
					}
					int[] temp = new int[acgt[k].length+1];
					for(int j=1;j<temp.length;j++){
						temp[j] = temp[j-1]+acgt[k][j-1];
					}
					acgt[k] = temp;
				}
				synchronized(lock3){
					lock3.decrement();
					lock3.notify();
				}
			}).start();
		}
		
		while(!lock3.finished()){
			System.err.println("waiting");
			synchronized(lock3){
				lock3.wait();
			}
		}
		
		
		System.err.println("prepare index: "+(System.currentTimeMillis()-time)/1000);
		
		
		time = System.currentTimeMillis();
		
		float[][] scores = new float[lists.length][];
		int numPer = scores.length/numThreads;

		Lock lock = new Lock(numThreads);
		QuickScanningSequenceScore[] lslims = new QuickScanningSequenceScore[numThreads];
		for(int i=0;i<lslims.length;i++){
			lslims[i] = (QuickScanningSequenceScore) lslim.clone();
		}
		
		
		
		
		//Thread[] threads = new Thread[numThreads];
		for(int i=0;i<numThreads;i++){
			final int j = i;
			new Thread(() -> {
				System.err.println("thread: "+j);
				for(int l=0;l<assign[j].length();l++){
					int k = assign[j].get(l);
					if((k)%(1000) == 0){
						System.err.println(k);
					}
					scores[k] = new float[lists[k].length()];
					try {
						computeScores(lslims[j],lists[k],k,prefK,scores[k]);
					} catch (WrongAlphabetException | WrongSequenceTypeException e) {
						e.printStackTrace();
						throw new RuntimeException(e);
					}
				}
				synchronized(lock){
					lock.decrement();
					lock.notify();
				}
			}).start();
		}
		

		while(!lock.finished()){
			System.err.println("waiting");
			synchronized(lock){
				lock.wait();
			}
		}
		/*System.out.println("test score: "+scores[0][10]);
		for(int i=0;i<scores.length;i++){
			if(scores[i] == null){
				throw new RuntimeException("null");
			}
		}*/
		
		/*for(int i=0;i<lists.length;i++){
			scores[i] = new float[lists[i].length()];
			computeScores(lslim,lists[i],i, prefK, scores[i]);
		}*/
		
		System.err.println("scores: "+(System.currentTimeMillis()-time)/1000);
		
		
		
		time = System.currentTimeMillis();
		
		read = new BufferedReader(new FileReader(genome));
		lastHeader = new StringBuffer();
		
		
		//numThreads = Math.min(2, numThreads);
		
		Lock lock2 = new Lock(numThreads);
		for(int i=0;i<numThreads;i++){
			lock2.setFinished(i, true);
		}

		
		BufferedOutputStream out = new BufferedOutputStream(new FileOutputStream(java.io.FileDescriptor.out),100000);
		
		LinkedList<ScoreRegion> scoreList = new LinkedList<>();

		int finished = 0;
		int totalSequenceIndex = -1;
		int lastPrinted = -1;
		while( (pair = LargeSequenceReader.readNextSequences(read, lastHeader, lslim.getLength())) != null ){

			

			ArrayList<Sequence> seqs = pair.getSecondElement();
			IntList starts = pair.getFirstElement();
			
				Iterator<Sequence> it = seqs.iterator();
				int itIdx = 0;
				
				while( it.hasNext() ) {
					final Sequence seq2 = it.next();
					
					totalSequenceIndex++;
					String id = seq2.getSequenceAnnotationByType("id", 0).getIdentifier().trim();
					int off = starts.get(itIdx);
					System.err.println(id);
					
					itIdx++;

					
					final int threadID = finished;
					final int myIndex = totalSequenceIndex;
					
					Thread temp = new Thread( ()->{
						Sequence seq = seq2;
						
						float[][] myScores = new float[2][seq.getLength()-lslim.getLength()+1];
						
						for(int d=0;d<2;d++){

							int idx1 = 0;
							for(int i=0;i<pow1.length-1;i++){
								idx1 += pow1[i+1]*seq.discreteVal(i);
							}
							int idx2 = 0;
							for(int i=0;i<pow2.length-1;i++){
								idx2 += pow2[i+1]*seq.discreteVal(i+prefK);
							}


							for(int j=0;j<seq.getLength()-lslim.getLength()+1;j++){
								idx1 = (idx1%pow1[0])*4 + seq.discreteVal(j+prefK-1);
								idx2 = (idx2%pow2[0])*4 + seq.discreteVal(j+prefK+restK-1);

								int sym = idx2/pow2[idxK-1];

								int idx3 =  lists[idx1].binarySearch(idx2, acgt[idx1][sym], acgt[idx1][sym+1]);//getIndex(lists[idx1], idx2, maxValue);

								float score = scores[idx1][idx3];
								
								myScores[d][ d==0 ? j : seq.getLength()-lslim.getLength()-j ] = score;
								
							}


							try {
								seq = seq.reverseComplement();
							} catch (OperationNotSupportedException e) {
								// TODO Auto-generated catch block
								e.printStackTrace();
								throw new RuntimeException();
							}
						}
						
						ScoreRegion sr = new ScoreRegion(myScores, myIndex, id, off);
						
						synchronized(lock2){
							add(scoreList,sr);
							lock2.setFinished(threadID, true);
							lock2.notify();
						}
					});

					synchronized(lock2){
						lock2.setFinished(finished, false);
						temp.start();
					}
					
					
					while(( finished = lock2.getFinished()) < 0){
						synchronized(lock2){
							lock2.wait();
						}
					}
					
					while(true){
						int lp2 = lastPrinted;
						lastPrinted = print(scoreList, lastPrinted, lock2,out);
						synchronized(lock2){
							if(scoreList.size() > numThreads && scoreList.getFirst().index != lastPrinted+1 ){

								lock2.wait();

							}else{
								break;
							}
						}
					}
				}

		}
		
		while(!lock2.allFinished()){
			synchronized(lock2){
				lock2.wait();
			}
		}
		
		lastPrinted = print(scoreList, lastPrinted, lock2,out);
		
		
		read.close();
		out.flush();
		out.close();
		
		System.err.println("scan: "+(System.currentTimeMillis()-time)/1000);

	}

	
	
	/*private static int getIndex2(IntList intList, int idx2, double maxValue){
		return intList.interpolationSearch(idx2, 0, intList.length());
	}



	private static int getIndex(IntList intList, int idx2, double maxValue) {
		
		int idx = (int) ( (idx2/maxValue) * intList.length() );
		
		int ill = intList.length();
		
		int start = 0;
		int end = ill;
		
		if(intList.get(idx) <= idx2){
			start = idx;
			if( intList.get((start+ill)/2)>idx2 ){
				end = (start+ill)/2;
			}
		}else{
			end = idx;
			if(intList.get(end/2)<=idx2){
				start = end/2;
			}
			
		}
		//System.out.println(start+" "+end);
		return intList.binarySearch(idx2, start, end);

	}
*/



	private boolean check(LinkedList<ScoreRegion> scoreList, int numThreads) {
		System.out.println("check");
		return scoreList.size() > numThreads*2;
	}

	private int print(LinkedList<ScoreRegion> scoreList, int lastPrinted, Lock lock, BufferedOutputStream out) throws IOException {
		int start = lastPrinted;
		
		LinkedList<ScoreRegion> regs2 = new LinkedList<>();
		
		synchronized(lock){
			Iterator<ScoreRegion> it = scoreList.iterator();
			while(it.hasNext()){
				ScoreRegion next = it.next();
				System.err.println(next.index+" <-> "+(lastPrinted+1));
				if(next.index == lastPrinted+1){
					regs2.add(next);
					lastPrinted++;
				}else{
					break;
				}
			}
			for(int i=start;i<lastPrinted;i++){
				scoreList.remove();
			}
		}
		
		
		Iterator<ScoreRegion> it = regs2.iterator();
		while(it.hasNext()){
			ScoreRegion next = it.next();
			StringBuilder sb = new StringBuilder();
			String tab = "\t";
			for(int i=0;i<next.scores[0].length;i++){
				sb.append(next.id);
				sb.append(tab);
				sb.append(next.offset+i);
				sb.append(tab);
				sb.append(next.scores[0][i]);
				sb.append(tab);
				sb.append(next.scores[1][i]);
				sb.append("\n");
				if((i+1)%100==0){
					out.write(sb.toString().getBytes());
					sb.delete(0, sb.length());
				}
			}
			out.write(sb.toString().getBytes());
			next.scores = null;
			out.flush();
		}
		regs2.clear();
		return lastPrinted;
	}

	private void add(LinkedList<ScoreRegion> scoreList, ScoreRegion sr) {
		Iterator<ScoreRegion> it = scoreList.iterator();
		int i=0;
		while(it.hasNext()){
			ScoreRegion temp = it.next();
			if(temp.index>sr.index){
				scoreList.add(i, sr);
				return;
			}else{
				i++;
			}
		}
		scoreList.add(i, sr);
	}

	private static IntList condense(IntList list) {
		
		int i=0, n=0;
		while(i<list.length()){
			while(i+1<list.length() && list.get(i)==list.get(i+1)){
				i++;
			}
			n++;
			i++;
		}
		
		IntList li = new IntList(n+1);
		i=0;
		while(i<list.length()){
			while(i+1<list.length() && list.get(i)==list.get(i+1)){
				i++;
			}
			li.add(list.get(i));
			i++;
		}
		return li;
	}




	private void computeScores(QuickScanningSequenceScore lslim, IntList list, int prefix, int prefixLength, float[] scoreAr) throws WrongAlphabetException, WrongSequenceTypeException {
		//System.out.println("prefix: "+prefix);
		int[] seq = new int[lslim.getLength()];
		
		DiscreteAlphabet al = (DiscreteAlphabet) lslim.getAlphabetContainer().getAlphabetAt(0);
		
		String[] syms = new String[(int) al.length()];
		for(int i=0;i<syms.length;i++){
			syms[i] = al.getSymbolAt(i);
		}
		
		double[] scores = new double[seq.length];
		
		fill(seq, prefix, prefixLength);
		
		int[] prevSeq = seq.clone();
		Arrays.fill(prevSeq, prefixLength,prevSeq.length,-1);
		
		lslim.fillInfixScore(seq,0,prefixLength, scores);
		
		
		int i=0;
		while(i<list.length()){
			
			int curr = list.get(i);
			
			//update(seq,curr,prefixLength);
			
			for(int k=prefixLength;k<seq.length;k++){
				seq[k] = curr/pow2[k-prefixLength];
				curr %= pow2[k-prefixLength];
			}
			
			//int offset = offset(prevSeq,seq,prefixLength);
			int offset=prefixLength;
			while(offset<seq.length && prevSeq[offset]==seq[offset]){
				offset++;
			}
			offset -= prefixLength;
			
			lslim.fillInfixScore(seq,prefixLength+offset,seq.length-prefixLength-offset, scores);
			
			double fullScore = ToolBox.sum(scores);
			
			scoreAr[i] = (float) fullScore;
			
			System.arraycopy(seq, prefixLength+offset, prevSeq, prefixLength+offset, seq.length-prefixLength-offset);
			
			i++;
		}
		
	}




	private void fill(int[] seq, int prefix, int len) {
		for(int i=0;i<len;i++){
			seq[i] = prefix/pow1[i];
			prefix %= pow1[i];
		}
	}



	
}
