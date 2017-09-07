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

package de.jstacs.clustering.distances;

import java.util.LinkedList;

import javax.naming.OperationNotSupportedException;

import de.jstacs.data.DeBruijnGraphSequenceGenerator;
import de.jstacs.data.WrongAlphabetException;
import de.jstacs.data.alphabets.DiscreteAlphabet;
import de.jstacs.data.sequences.CyclicSequenceAdaptor;
import de.jstacs.data.sequences.WrongSequenceTypeException;
import de.jstacs.sequenceScores.statisticalModels.StatisticalModel;
import de.jstacs.utils.Pair;
import de.jstacs.utils.RealTime;

/**
 * Class for a distance metric between {@link StatisticalModel}s based on the correlation of score
 * profiles on De Bruijn sequences.
 * 
 * For two {@link StatisticalModel}s {@latex.inline $M_1$} and {@latex.inline $M_2$}, we compute the score profiles
 * {@latex.inline $s_1(x,M_1)$} and {@latex.inline $s_2(x,M_2)$} on a De Bruijn sequence {@latex.inline x} of length {@latex.inline $|A|^n$}. 
 * The distance is then defined based on the Pearson correlation as {@latex.inline $1 - cor( s_1(x,M_1), s_2(x,M_2) )$} between these score profiles,
 * maximizing over suitable shifts of the score profiles and both strand orientations. 
 * 
 * @author Jan Grau
 *
 */
public class SequenceScoreDistance extends DistanceMetric<StatisticalModel> {

	/**
	 * The De Bruijn sequences
	 */
	protected CyclicSequenceAdaptor[] seqs;
	
	/**
	 * if exponential scores should be used
	 */
	protected boolean exp;
	
	/**
	 * Creates a new distance.
	 * @param alphabet the alphabet of the models that may be compared
	 * @param n the length of n-mers represented in the De Bruijn sequence
	 * @param exp if exponential scores should be used
	 * @throws WrongAlphabetException if the sequence of this alphabet could not be created
	 * @throws WrongSequenceTypeException if the sequence could not be created
	 */
	public SequenceScoreDistance(DiscreteAlphabet alphabet, int n, boolean exp) throws WrongAlphabetException, WrongSequenceTypeException{
		this(DeBruijnGraphSequenceGenerator.generate( alphabet, n ), exp);
	}
	
	/**
	 * Creates a new distance for a given set of sequences.
	 * @param seqs the sequences
	 * @param exp if exponential scores should be used
	 */
	protected SequenceScoreDistance(CyclicSequenceAdaptor[] seqs, boolean exp){
		this.seqs = seqs;
		//System.out.println(seqs[0]);
		this.exp = exp;
	}
	
	/**
	 * Returns the score profile for the model.
	 * @param o the mode
	 * @param rc if the reverse complement should be considered
	 * @return the score profile
	 * @throws Exception if the score could not be computed
	 */
	public double[][] getProfile(StatisticalModel o, boolean rc) throws Exception{
		return DeBruijnMotifComparison.getProfilesForMotif( seqs, o, rc, exp );
	}
	
	/**
	 * Returns the distance between the two score profiles. 
	 * @param profiles1 the first profile
	 * @param profiles1Rc the reverse complementary version of the first profile
	 * @param profiles2 the second profile
	 * @param maxShift the maximum allowed shift between the profiles
	 * @return the distance
	 * @throws Exception if the distance could not be computed
	 */
	public double getDistance(double[][] profiles1, double[][] profiles1Rc, double[][] profiles2, int maxShift) throws Exception{
		Pair<Integer,Double> newFwd = DeBruijnMotifComparison.compare( profiles1[0], profiles2[0], maxShift ); 
		Pair<Integer,Double> newRc = DeBruijnMotifComparison.compare( profiles1Rc[0], profiles2[0], maxShift );
		
		return 1.0-Math.max( newFwd.getSecondElement(), newRc.getSecondElement() );
	}
	
	/**
	 * Returns the distance between a score profile and a model. 
	 * @param profiles1 the first profile
	 * @param profiles1Rc the reverse complementary version of the first profile
	 * @param o2 the model
	 * @param motif1Length the length of the motif used to compute the first profile
	 * @return the distance
	 * @throws Exception if the distance could not be computed
	 */
	public double getDistance(double[][] profiles1, double[][] profiles1Rc, StatisticalModel o2, int motif1Length) throws Exception{
		
		//int maxShift = Math.max( motif1Length - (int)Math.floor(o2.getLength()/3), o2.getLength() - (int)Math.floor( motif1Length/3 ) );
		int maxShift = (int)Math.ceil(Math.min(o2.getLength(), motif1Length)*2.0/3.0);
		double[][] profiles2 = getProfile( o2, false );
		//System.out.println(profiles1[0].length+" "+profiles2[0].length);
		return getDistance( profiles1, profiles1Rc, profiles2, maxShift );
	}
	
	@Override
	public double getDistance( StatisticalModel o1, StatisticalModel o2 ) throws Exception {
		
		int motif1Length = o1.getLength();
		
		double[][] profiles1 = getProfile( o1, false );
		double[][] profiles1Rc = getProfile( o1, true );
		
		return getDistance( profiles1, profiles1Rc, o2, motif1Length );
		
	}
	
	/**
	 * Multi-threaded computation of the pairwise distance matrix.
	 * @param numThreads the number of threads
	 * @param objects the models
	 * @return the distance matrix
	 * @throws Exception if the distance could not be computed
	 */
	public double[][] getPairwiseDistanceMatrix(int numThreads, StatisticalModel... objects) throws Exception {
		if(numThreads == 1){
			double[][] matrix = new double[objects.length][];
			for(int i=0;i<matrix.length;i++){
				matrix[i] = new double[i];
				double[][] prof1 = this.getProfile( objects[i], false );
				double[][] prof1Rc = this.getProfile( objects[i], true );
				fillDistanceRow( i, matrix[i], prof1, prof1Rc, this, objects );
			}
			return matrix;
		}else{
			LinkedList<Worker> available = new LinkedList<Worker>();
			Worker[] workers = new Worker[numThreads-1];
			for(int i=0;i<workers.length;i++){
				workers[i] = new Worker( objects, available, i );
				available.add( workers[i] );
				(new Thread(workers[i])).start();
			}
			RealTime rt = new RealTime();
			double[][] matrix = new double[objects.length][];
			for(int i=0;i<matrix.length;i++){
				System.out.println("main in row "+i+" "+rt.getElapsedTime());
				matrix[i] = new double[i];
				double[][] prof1 = null;
				double[][] prof1Rc = null;
				int length = 0;
				synchronized(objects[i]){
					prof1 = this.getProfile( objects[i], false );
					length = objects[i].getLength();
					prof1Rc = this.getProfile( objects[i], true );
				}
				int a = 0;
				while(a == 0){
					synchronized( available ) {
						a = available.size();
					}
					if(a == 0){
						synchronized(this){
					//		System.out.println("main waits");
							wait(1000);
						}
					}
				}
				Worker curr = null;
				synchronized( available ) {
					curr = available.pop();
				}
			//	System.out.println("main found worker "+curr.index);
				curr.setRow( matrix[i], prof1, prof1Rc, length );
			//	System.out.println("main set row in worker "+curr.index);
			}
			int a = 0;
			while(a < workers.length){
				synchronized( available ) {
					a = available.size();
				}
				if(a == workers.length){
					synchronized(this){
		//				System.out.println("main waits");
						wait(100);
					}
				}
			}
			
			for(int i=0;i<workers.length;i++){
				workers[i].stop();
			}
			return matrix;
		}
	}
	
	private static void fillDistanceRow(int i, double[] row, double[][] prof1, double[][] prof1Rc, SequenceScoreDistance metric, StatisticalModel... objects) throws Exception{
		for(int j=0;j<i;j++){
			row[j] = metric.getDistance( prof1, prof1Rc, objects[j], objects[i].getLength() );
		}
	}

	
	private class Worker implements Runnable{

		private int index;
		private boolean stop;
		private double[] row;
		private double[][] prof1;
		private double[][] prof1Rc;
		private SequenceScoreDistance metric;
		private StatisticalModel[] objects;
		private LinkedList<Worker> available;
		private int length;
		
		private Worker(StatisticalModel[] objects, LinkedList<Worker> available, int index){
			this.objects = objects;
			this.metric = SequenceScoreDistance.this;
			this.available = available;
			this.index = index;
		}
		
		private void setRow(double[] row, double[][] prof1, double[][] prof1Rc, int length){
			//System.out.println("set row");
			synchronized( this ) {
				this.row = row;
				this.prof1 = prof1;
				this.prof1Rc = prof1Rc;
				this.length = length;
				
				//System.out.println("notified this in set");
				notify();
			}
		}
		
		private void stop(){
			this.stop = true;
			synchronized(this){
				//System.out.println("notified this in stop");
				this.notify();
			}
		}
		
		@Override
		public void run() {
			synchronized(this){
				try {
					//System.out.println("worker waits");
					this.wait();
				} catch ( InterruptedException e ) {
					e.printStackTrace();
				}
			}
			while(!stop){
				
				if(row != null){
					//System.out.println("worker "+index+" works");
					try {
						fillDistanceRow( );
						synchronized(this){
							this.row = null;
							synchronized( available ) {
					//			System.out.println("worker "+index+" available");
								available.add( this );
							}
							synchronized( metric ) {
					//			System.out.println("worker "+index+" notified metric");
								metric.notify();
							}
							if(stop){
								return;
							}
							try {
					//			System.out.println("worker "+index+" waits");
								wait();
							} catch ( InterruptedException e ) {
								e.printStackTrace();
							}
						}
					} catch ( Exception e ) {
						e.printStackTrace();
						this.stop = true;
						throw new RuntimeException();
					}
				}
			}
		}
		
		private void fillDistanceRow() throws Exception{
			for(int j=0;j<row.length;j++){
				synchronized(objects[j]){
					row[j] = metric.getDistance( prof1, prof1Rc, objects[j], length );
				}
			}
		}
		
	}
}
