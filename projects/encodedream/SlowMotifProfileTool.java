package projects.encodedream;
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
import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.LinkedList;

import javax.naming.OperationNotSupportedException;

import de.jstacs.data.WrongAlphabetException;
import de.jstacs.data.alphabets.DiscreteAlphabet;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.data.sequences.WrongSequenceTypeException;
import de.jstacs.sequenceScores.QuickScanningSequenceScore;
import de.jstacs.utils.IntList;
import de.jstacs.utils.LargeSequenceReader;
import de.jstacs.utils.Pair;
import de.jstacs.utils.ToolBox;

public class SlowMotifProfileTool {

	
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
	
	
	
	public void run(QuickScanningSequenceScore lslim, String genome, int numThreads, BufferedOutputStream out) throws Exception {
		
		
		
		BufferedReader read = new BufferedReader(new FileReader(genome));
		StringBuffer lastHeader = new StringBuffer();
		
		long time = System.currentTimeMillis();

		Pair<IntList,ArrayList<Sequence>> pair = null;
		
		LinkedList<ScoreRegion> scoreList = new LinkedList<>();
		Lock lock2 = new Lock(numThreads);
		for(int i=0;i<numThreads;i++){
			lock2.setFinished(i, true);
		}

		int finished = 0;
		int totalSequenceIndex = -1;
		int lastPrinted = -1;
		while( (pair = LargeSequenceReader.readNextSequences(read, lastHeader, lslim.getLength(), lslim.getAlphabetContainer())) != null ){

			

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

							for(int j=0;j<seq.getLength()-lslim.getLength()+1;j++){

								float score = (float) lslim.getLogScoreFor(seq, j);
								
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
		//out.close();
		
		System.err.println("scan: "+(System.currentTimeMillis()-time)/1000);

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

	
}
