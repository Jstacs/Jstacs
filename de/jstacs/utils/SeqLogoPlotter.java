package de.jstacs.utils;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Font;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.Polygon;
import java.awt.RenderingHints;
import java.awt.Shape;
import java.awt.geom.AffineTransform;
import java.awt.geom.Area;
import java.awt.geom.Ellipse2D;
import java.awt.geom.Rectangle2D;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.LinkedList;

import javax.imageio.ImageIO;

import cern.jet.stat.Gamma;
import de.jstacs.algorithms.optimization.termination.SmallDifferenceOfFunctionEvaluationsCondition;
import de.jstacs.data.DataSet;
import de.jstacs.data.EmptyDataSetException;
import de.jstacs.data.WrongAlphabetException;
import de.jstacs.data.alphabets.DNAAlphabetContainer;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.io.ArrayHandler;
import de.jstacs.io.NonParsableException;
import de.jstacs.io.XMLParser;
import de.jstacs.results.PlotGeneratorResult.PlotGenerator;
import de.jstacs.sequenceScores.statisticalModels.differentiable.directedGraphicalModels.structureLearning.measures.btMeasures.BTExplainingAwayResidual;
import de.jstacs.sequenceScores.statisticalModels.trainable.TrainableStatisticalModel;
import de.jstacs.sequenceScores.statisticalModels.trainable.TrainableStatisticalModelFactory;
import de.jstacs.sequenceScores.statisticalModels.trainable.mixture.AbstractMixtureTrainSM.Parameterization;
import de.jstacs.sequenceScores.statisticalModels.trainable.mixture.MixtureTrainSM;
import de.jstacs.utils.graphics.GraphicsAdaptor;

/**
 * Class with static methods for plotting sequence logos of DNA motifs, i.e., position weight matrices defined over a {@link DNAAlphabet}.
 * In general, sequence logos can be plotted to any {@link Graphics2D} object, e.g., for on-screen printing using the <code>plotLogo</code> methods.
 * 
 * For convenience, the method {@link SeqLogoPlotter2#plotLogoToPNG(String, int, double[][])} can be used to directly plot a sequence logo
 * to a PNG file with a given height and automatically chosen aspect ratio.
 * 
 * @author Jan Grau
 *
 */
public class SeqLogoPlotter {

	/**
	 * {@link PlotGenerator} for plotting sequence logos.
	 * 
	 * @author Jan Grau
	 *
	 */
	public static class SeqLogoPlotGenerator implements PlotGenerator{

		private double[][] pwm;
		private int height;
		
		/**
		 * Creates a new {@link SeqLogoPlotGenerator} for the given PWM using the specified height of the plot.
		 * The corresponding width is computed automatically using the {@link SeqLogoPlotter#getWidth(int, double[][])} method.
		 * 
		 * @param pwm the PWM
		 * @param height the height of the plot
		 */
		public SeqLogoPlotGenerator(double[][] pwm, int height){
			this.pwm = pwm;
			this.height = height;
		}
		
		/**
		 * Creates a {@link SeqLogoPlotGenerator} from its XML representation.
		 * @param xml the XML representation
		 * @throws NonParsableException if XML could not be parsed
		 */
		public SeqLogoPlotGenerator(StringBuffer xml) throws NonParsableException{
			this.pwm = (double[][]) XMLParser.extractObjectForTags(xml, "pwm");
			this.height = (Integer) XMLParser.extractObjectForTags(xml, "height");
		}
		
		@Override
		public StringBuffer toXML() {
			StringBuffer xml = new StringBuffer();
			XMLParser.appendObjectWithTags(xml, pwm, "pwm");
			XMLParser.appendObjectWithTags(xml, height, "height");
			return xml;
		}

		@Override
		public void generatePlot(GraphicsAdaptor ga) throws Exception {
			
			int width = SeqLogoPlotter.getWidth(height, pwm);
			
			SeqLogoPlotter.plotLogo(ga.getGraphics(width, height), height, pwm);
			
		}
		
	}

	
	private static Color getColor(double[] pwm, float ic){
		float[] a = new float[]{0f,1f,0f};
		float[] c = new float[]{0f,0f,1f};
		//float[] g = new float[]{1f,1f,0f};
		float[] g = Color.ORANGE.getColorComponents( null );
		float[] t = new float[]{1f,0f,0f};
		
		
		
		float[] nc = new float[3];
		for(int i=0;i<3;i++){
			nc[i] = (float)( pwm[0]*a[i] + pwm[1]*c[i] + pwm[2]*g[i] + pwm[3]*t[i] );
		}
		return new Color( nc[0], nc[1], nc[2], ic );
	}
	
	public static double[][] getDeps(Sequence[] seqs, double[] w){
	
		double[][][][] counts = new double[seqs[0].getLength()][seqs[0].getLength()][4][4];
		
		
		double sum = 0;
		
		for(int i=0;i<seqs.length;i++){
			double weight = (w == null ? 1d : w[i]);
			sum += weight;
			for(int k=0;k<seqs[i].getLength();k++){
				//counts01[k][seqs[i].discreteVal( k )] += weight;
				for(int l=0;l<seqs[i].getLength();l++){
					counts[k][l][seqs[i].discreteVal( k )][seqs[i].discreteVal( l )]+=weight;
					
				}
			}
			
			
			
		}
		
		double[][] mis = new double[seqs[0].getLength()][seqs[0].getLength()];
		
		
		for(int k=1;k<seqs[0].getLength();k++){
			for(int l=0;l<k;l++){
				//if(exclude == null || (!exclude[k] && !exclude[l]) ){		
					double mi = 0;
					for(int a=0;a<4;a++){
						for(int b=0;b<4;b++){
							if(counts[k][l][a][b]>0){
								mi += counts[k][l][a][b] * ( Math.log( counts[k][l][a][b] ) - Math.log(counts[k][k][a][a]) - Math.log( counts[l][l][b][b] ) + Math.log( sum ) );
							}
						}
					}
					mis[k][l] = mi;
					mis[l][k] = mi;
				//}
			}
		}
		
		return mis;
		
	}
	
	public static int getHeightForColorLogo(int numSeqs, int numOne, int numPerChunk, int oneHeight, int blockSpacer){
		
		return (numSeqs/numOne*oneHeight) + (numSeqs/numPerChunk*blockSpacer);
		
	}
	
	public static int getHeightForDependencyLogo( int seqLength, int numSeqs, int[] chunkHeights, int width, int blockSpacer){
		
		int height = getHeightForColorLogo( numSeqs, chunkHeights, blockSpacer );
		
		int topMargin = (int)width/5;
		
		height += 1.5*topMargin;
		
		return height;
	}
	
	public static int getHeightForColorLogo(int numSeqs, int[] chunkHeights, int blockSpacer){
		
		int height = 0;
		
		for(int i=0;i<chunkHeights.length;i++){
			height += chunkHeights[i] + blockSpacer;
		}
		return height;
	}
	
	private static void plotScale(Graphics g, int offx, int offy, int height){
		g = g.create();
		g.setColor( Color.BLACK );
		Rectangle2D rect = g.getFontMetrics().getStringBounds( "2", g ); 
		int w = (int)rect.getWidth();
		int h = -(int)rect.getCenterY();
		g.drawLine( offx, offy, offx, offy+height );
		g.drawString( "2", offx-2*w, offy+h );
		g.drawLine( offx-(int)(0.8*w), offy, offx, offy );
		g.drawString( "1", offx-2*w, (int)(offy+0.5*height)+h );
		g.drawLine( offx-(int)(0.8*w), offy+(int)(0.5*height), offx, offy+(int)(0.5*height) );
		g.drawString( "0", offx-2*w, (int)(offy+height)+h );
		g.drawLine( offx-(int)(0.8*w), offy+height, offx, offy+height );
	}
	
	
	public static int[] getBorders(DataSet data, int[][] minmax, int[] steps) throws Exception{
		
		int off = 0;
		int[] maxBords = new int[minmax.length];
		for(int i=0;i<minmax.length;i++){
			double max = Double.NEGATIVE_INFINITY;
			int maxBord = -1;
			for(int j=minmax[i][0];j<=minmax[i][1] && off+j < data.getNumberOfElements();j+=steps[i]){

				Sequence[] sub1 = new Sequence[j];
				System.arraycopy( data.getAllElements(), off, sub1, 0, sub1.length );
				Sequence[] sub2 = new Sequence[data.getNumberOfElements()-(off+j)];
				System.arraycopy( data.getAllElements(), off+j, sub2, 0, sub2.length );
				TrainableStatisticalModel pwm1 = TrainableStatisticalModelFactory.createPWM( DNAAlphabetContainer.SINGLETON, data.getElementLength(), 0.0 );
				pwm1.train( new DataSet( "", sub1 ) );
				TrainableStatisticalModel pwm2 = TrainableStatisticalModelFactory.createPWM( DNAAlphabetContainer.SINGLETON, data.getElementLength(), 0.0 );
				pwm2.train( new DataSet("", sub2) );
				double[] weights = new double[]{sub1.length/(double)(sub1.length+sub2.length),sub2.length/(double)(sub1.length+sub2.length)};
				//double[] weights = new double[]{0.5,0.5};
				//System.out.println(Arrays.toString( weights ));
				MixtureTrainSM mtsm = new MixtureTrainSM( pwm1.getLength(), new TrainableStatisticalModel[]{pwm1,pwm2},weights, 1, 1, new SmallDifferenceOfFunctionEvaluationsCondition( 1E-6 ), Parameterization.THETA);
				double ll = 0;
				for(int k=off;k<data.getNumberOfElements();k++){
					ll += mtsm.getLogProbFor( data.getElementAt( k ) );

				}
				if(ll > max){
					max = ll;
					maxBord = j;
				}

			}
			maxBords[i] = maxBord;
			off += maxBord;
		}
		return maxBords;
	}
	
	public static BufferedImage plotDefaultDependencyLogoToBufferedImage(DataSet data, double[] weights, int width) throws Exception{
		
		
		int[] numPerChunk = new int[]{Math.min( 250, (int)Math.round( data.getNumberOfElements()*0.1 ) ), Math.min( 1250, (int)Math.round( data.getNumberOfElements()*0.3 ) ),0};
		numPerChunk[2] = data.getNumberOfElements() - numPerChunk[0] - numPerChunk[1];
		
		int logoHeight = (int)Math.round( width/(double)10 );
		
		int[] chunkHeights = new int[]{60,75,150};
		
		int height = getHeightForDependencyLogo( data.getElementLength(), data.getNumberOfElements(), chunkHeights, width, logoHeight );//getHeightForColorLogo( data.getNumberOfElements(), numOne, numPerChunk, oneHeight, logoHeight );
		
		BufferedImage img = new BufferedImage( width, height, BufferedImage.TYPE_INT_RGB);
		Graphics2D graph = (Graphics2D)img.getGraphics();
		
		graph.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
		graph.setRenderingHint(RenderingHints.KEY_RENDERING, RenderingHints.VALUE_RENDER_QUALITY);
		
		plotDependencyLogo( data, null, 1, null, weights, graph, width, 0, 0, numPerChunk,chunkHeights, 0.03, logoHeight, false, 3, false, true, true, 0.1 );
		
		return img;
	}
	
	public static int plotDependencyLogo(DataSet seqs, Object[] labels, int ticPeriod, double[][] classProbs, double[] weights, Graphics2D graph, int width, int offx, int offy, int[] numPerChunk, int[] chunkHeights, double minPercent, int logoHeight, boolean highlightMaxDeps, int numBestForSorting, boolean sortGlobally, boolean sortByWeights, boolean scaleByDeps, double threshold) throws Exception {

		int totalWidth = width;
		
		graph = (Graphics2D)graph.create();
		
		int leftMargin = width/15;
		
		int tempWidth = width - 2*leftMargin;
		
		if(labels == null){
			labels = new String[seqs.getElementLength()];
			for(int i=0;i<labels.length;i++){
				labels[i] = (i+1)+"";
			}
		}
		int partWidth = (int)Math.floor( (double)tempWidth/(double)(labels.length) );
		
		
		Font font = new Font(graph.getFont().getName(),Font.BOLD,width/30);
		graph.setFont( font );
		
		int fac = 30;
		
		double maxWidth = 0;
		for(int i=0;i<labels.length;i++){
			double temp = graph.getFontMetrics().getStringBounds( labels[i].toString(), graph ).getWidth();
			if(temp > maxWidth){
				maxWidth = temp;
			}
		}
		
		if( maxWidth > partWidth*0.9*ticPeriod ){
			fac *= maxWidth/(partWidth*0.9*ticPeriod);
			font = new Font(graph.getFont().getName(),Font.BOLD,width/fac);
			graph.setFont( font );
		}
		
		graph.setStroke( new BasicStroke( width/300f ) );
		
		
		
		double symHeight = graph.getFontMetrics().getStringBounds( labels[labels.length-1].toString(), graph ).getHeight();
		
		
		//int topMargin = (int)(symHeight*6);
		int topMargin = (int)width/5;
		
		graph.setColor( Color.WHITE );
		graph.fillRect( offx, offy, width, getHeightForColorLogo( seqs.getNumberOfElements(), chunkHeights, logoHeight ) + topMargin*2 );

		
		double labWidth = graph.getFontMetrics().getStringBounds( "0.000E000", graph ).getWidth();
		
		
		width -= 2*leftMargin;
		
		

		
		
		
		Sequence[] seqs2 = seqs.getAllElements();
		
		int[] heights = plotColorLogo( seqs2, weights, graph, width, offx+leftMargin, offy+topMargin, numPerChunk, chunkHeights, minPercent, logoHeight, numBestForSorting, sortGlobally, sortByWeights, threshold );
		
		seqs = new DataSet("",seqs2);
		
		graph.setColor( Color.BLACK );
		int off = topMargin;
		for(int i=0;i<heights.length;i++){
			
			String lab = numPerChunk[i]+"";
			

			Rectangle2D rect = graph.getFontMetrics().getStringBounds( lab, graph );
			AffineTransform back = graph.getTransform();
			graph.rotate( -Math.PI/2 );
			//System.out.println("LAB: "+lab+" "+(-(int)((heights[i]+off-logoHeight)/2 + rect.getCenterX()))+" "+((int)(leftMargin-rect.getHeight()/2d)));
			graph.drawString(lab,-(int)((heights[i]+off-logoHeight)/2 + rect.getCenterX()), (int)(leftMargin-rect.getHeight()/2d));
			graph.setTransform( back );
			
			off = heights[i];
		}
		
		
		off = partWidth/2;
		
		addLabels(labels,ticPeriod,graph,partWidth,offx,leftMargin,off,symHeight,topMargin,offy,heights[ heights.length-1 ]);

		double[][] deps = null;
		if(classProbs == null || classProbs.length==1){
			deps = getDeps( seqs2, classProbs == null ? null : classProbs[0] );
		}else{
			BTExplainingAwayResidual ear = new BTExplainingAwayResidual( new double[]{0,0} );
			deps = ear.getEAR( new DataSet("",seqs2), new DataSet("",seqs2), classProbs[0], classProbs[1], seqs2[0].getLength() );
		}
		boolean[][] isMax = new boolean[deps.length][deps.length];

		
		if(highlightMaxDeps){
		
			for(int i=0;i<deps.length;i++){
				double max = ToolBox.max( deps[i] );
				for(int j=0;j<deps[i].length;j++){
					if(deps[i][j] == max){
						isMax[i][j] = true;
					}
				}
			}
			
			for(int i=1;i<deps.length;i++){
				for(int j=0;j<i;j++){
					isMax[i][j] = isMax[i][j] || isMax[j][i];
					isMax[j][i] = isMax[i][j];
				}
			}
			
		}
		
		
		double[][] ps = new double[deps.length][deps[0].length];
		double f = deps.length*(deps.length-1)/2.0;
		double maxp = -Math.log10( 1E-300/f ), minp = -Math.log10( 0.01/f );
		for(int i=1;i<deps.length;i++){
			for(int j=0;j<i;j++){
				double p = -Math.log10( Gamma.incompleteGammaComplement( 9.0/2.0, deps[i][j]*2.0/2.0 ) );
				if(p > maxp){
					p = maxp;
				}
				ps[i][j] = p;
			/*	if(p > maxp){
					maxp = p;
				}
				/*if(p < minp){
					minp = p;
				}*/
			}
		}
		
		if(scaleByDeps){
			
			for(int i=1;i<deps.length;i++){
				for(int j=0;j<i;j++){
					if(ps[i][j]>minp){
						ps[i][j] = Math.log(deps[i][j]/( classProbs == null ? (double)seqs2.length : ToolBox.sum( classProbs[0] )));
					}else{
						ps[i][j] = Double.NEGATIVE_INFINITY;
					}
					/*if(ps[i][j] > maxp){
						maxp = ps[i][j];
					}*/
				}
			}
			minp = Math.log(0.01);//Math.sqrt( threshold );
			maxp = Math.log( Math.log( 4.0 ) );//Math.log( 4.0 )/2.0;
		}
		
		
		//System.out.println("max: "+maxp+", min: "+minp);
		
		
		graph.setColor( Color.BLACK );
		for(int i=1;i<deps.length;i++){
			for(int j=0;j<i;j++){
				if(ps[i][j] > minp && !isMax[i][j]){
					Color c = new Color( 0f, 0f, 0f, (float)((ps[i][j] - minp)/(maxp-minp)) );
					//System.out.println(c.getAlpha());
					graph.setColor( c );
					int h = (int)((topMargin-symHeight*1.5)*(i-j)/(double)deps.length);
					graph.drawArc( offx+leftMargin+off+j*partWidth,offy+(int)(topMargin-symHeight*1.5)-h, (i-j)*partWidth, h*2, 0, 180 );
				}
			}
		}
		for(int i=1;i<deps.length;i++){
			for(int j=0;j<i;j++){
				if(ps[i][j] > minp && isMax[i][j]){
					Color c = new Color( 1f, 0f, 0f, (float)((ps[i][j] - minp)/(maxp-minp)) );
					//System.out.println(c.getAlpha());
					graph.setColor( c );
					int h = (int)((topMargin-symHeight*1.5)*(i-j)/(double)deps.length);
					graph.drawArc( offx+leftMargin+off+j*partWidth,offy+(int)(topMargin-symHeight*1.5)-h, (i-j)*partWidth, h*2, 0, 180 );
				}
			}
		}
		
		
		
		off = 0;
		for(int i=0;i<heights.length;i++){
			
			double[][] pwm = PFMComparator.getPWM( seqs, off, off + numPerChunk[i] );
			off += numPerChunk[i];
			
			for(int j=0;j<pwm.length;j++){
				plotLogo( graph, offx+leftMargin+j*partWidth, heights[i], partWidth, logoHeight, pwm[j] );
			}
			plotScale( graph, offx+leftMargin, heights[i]-logoHeight, logoHeight-1 );
			
		}
		
		
		return heights[heights.length-1];
	}
	
	public static void addLabels(Object[] labels, int ticPeriod, Graphics2D graph, int partWidth, int offx, int leftMargin, int off, double symHeight, int topMargin,int offy, int lastHeight) {
		graph.setColor( Color.BLACK );
		for(int i=0;i<labels.length;i++){
			int per = (i+ticPeriod-1)%ticPeriod;
			int c = (int)graph.getFontMetrics().getStringBounds( labels[i].toString(), graph ).getCenterX();
			if(labels[i] instanceof Integer){
				int num = (Integer) labels[i];
				if(num < 0){
					c += (int)graph.getFontMetrics().getStringBounds( "-", graph ).getCenterX();
				}
				per = ( ( num < 0 ? -num : num )+ticPeriod)%ticPeriod;

			}
			if(  per == 0 ){
				if(ticPeriod != 1){
					graph.drawString( labels[i].toString(), offx+leftMargin+off+i*partWidth - c, offy+(int)(topMargin-symHeight/1.5) );
				}else{
					graph.drawString( labels[i].toString(), offx+leftMargin+off+i*partWidth - c, offy+(int)(topMargin-symHeight/2) );
				}
				graph.drawString( labels[i].toString(), offx+leftMargin+off+i*partWidth - c, lastHeight+(int)(symHeight*2.5) );
				if(ticPeriod != 1){
					graph.drawLine( offx+leftMargin+off+i*partWidth, offy+(int)(topMargin-symHeight/2), offx+leftMargin+off+i*partWidth, offy+(int)(topMargin-symHeight/2+symHeight*0.4) );
				}
				graph.drawLine( offx+leftMargin+off+i*partWidth, lastHeight+(int)symHeight, offx+leftMargin+off+i*partWidth, lastHeight+(int)(1.5*symHeight) );
			}else{
				if(ticPeriod != 1){
					graph.drawLine( offx+leftMargin+off+i*partWidth, offy+(int)(topMargin-symHeight/2+symHeight*0.2), offx+leftMargin+off+i*partWidth, offy+(int)(topMargin-symHeight/2+symHeight*0.4) );
				}

				graph.drawLine( offx+leftMargin+off+i*partWidth, lastHeight+(int)symHeight, offx+leftMargin+off+i*partWidth, lastHeight+(int)(1.25*symHeight) );
			}
		}

		graph.drawLine( offx+leftMargin+off, lastHeight+(int)symHeight, offx+leftMargin+off+(labels.length-1)*partWidth, lastHeight+(int)symHeight );

		if(ticPeriod != 1){
			graph.drawLine( offx+leftMargin+off, offy+(int)(topMargin-symHeight/2+symHeight*0.4), offx+leftMargin+off+(labels.length-1)*partWidth, offy+(int)(topMargin-symHeight/2+symHeight*0.4) );
		}
		
	}

	public static int plotColorLogo(Sequence[] seqs, double[] weights, Graphics2D graph, int width, int offx, int offy, int numOne, int numPerChunk, int oneHeight, int blockSpacer, int numBestForSorting, boolean sortGlobally, boolean sortByWeights, double threshold) throws Exception{
		int[] nums = new int[(int)Math.ceil( seqs.length/(double)numPerChunk )];
		int[] oneNums = new int[nums.length];
		for(int i=0,k=0;i<seqs.length;i+=numPerChunk,k++){
			if(i + numPerChunk <= seqs.length){
				nums[k] = numPerChunk;
			}else{
				int numCurr = (int) Math.floor( (seqs.length-i)/numOne ) * numOne;
				nums[k] = numCurr;
			}
			oneNums[k] = numOne;
		}
		//System.out.println(Arrays.toString( nums ));
		int[] offs = plotColorLogo( seqs, weights, graph, width, offx, offy, oneNums, nums, oneHeight, blockSpacer, numBestForSorting, sortGlobally, sortByWeights, threshold );
		return offs[offs.length-1];
	}
	
	public static int[] plotColorLogo(Sequence[] seqs, double[] weights, Graphics2D graph, int width, int offx, int offy, int[] numPerChunk, int[] chunkHeights, double minPercent, int blockSpacer, int numBestForSorting, boolean sortGlobally, boolean sortByWeights, double threshold) throws Exception{
		
		int[] offs = new int[numPerChunk.length];
		
		//double[] tempWeights = null;
		
		if(weights == null){
			weights = new double[seqs.length];
			Arrays.fill( weights, 1.0 );
		}else{

			/*double max = ToolBox.max( weights );
			double min = ToolBox.min( weights );
			for(int i=0;i<weights.length;i++){
				weights[i] = (weights[i] - min)/(max-min);
			}*/
			ComparableElement<Sequence, Double>[] els = new ComparableElement[seqs.length];
			for(int i=0;i<els.length;i++){
				els[i] = new ComparableElement<Sequence, Double>( seqs[i], weights[i] );
			}
			Arrays.sort( els );
			
			for(int i=0;i<els.length;i++){
				seqs[i] = els[els.length-1-i].getElement();
				weights[i] = els[els.length-1-i].getWeight();
			}
		}
		
		int off = 0;
		for(int i=0;i<numPerChunk.length;i++){
			//System.out.println("PLOTTING for "+i);
			SeqLogoPlotter.plotColorLogo( graph, seqs, weights, off, off+numPerChunk[i], width, offx, offy, chunkHeights[i], minPercent, numBestForSorting, sortGlobally, sortByWeights, threshold );
			offy += chunkHeights[i] + blockSpacer;
			Color back = graph.getColor();
			graph.setColor( Color.WHITE );
			graph.fillRect( offx, offy-blockSpacer, width, blockSpacer );
			graph.setColor( back );
			off += numPerChunk[i];
			offs[i] = offy;
		}
		
		return offs;
	}
	
	
	public static int plotColorLogo(Graphics2D graph, Sequence[] seqs, double[] weights, int start, int end, int width, int offx, int offy, int chunkHeight, double minPercent, int numBestForSorting, boolean sortGlobally, boolean sortByWeights, double threshold) throws Exception {
		
		double mi = ToolBox.min( weights );
		double ma = ToolBox.max( weights );
		mi -= (ma-mi)*1E-6;
		if(ma==mi){
			mi=0;
		}
		
		graph = (Graphics2D)graph.create();
		
		Sequence[] temp = new Sequence[end-start];
		System.arraycopy( seqs, start, temp, 0, end-start );
		
		Pair<Sequence,Double>[] sortTemp = new Pair[end-start];
		for(int i=0;i<sortTemp.length;i++){
			sortTemp[i] = new Pair<Sequence, Double>( seqs[start+i], weights[start+i] );
		}
		
		
		//double[] vals = getInformation( temp, numBestForSorting );
		
		//final int[] ord = ToolBox.order( vals, true );
		


		double[][] fullPFM = PFMComparator.getPFM( new DataSet("",temp) );

		int numSort = 6;
		
		
		LinkedList<ComparableElement<Integer, Double>> pivots = new LinkedList<ComparableElement<Integer,Double>>();

		Pair<Sequence,Double>[][] sorted = sortLocal2(sortTemp, numBestForSorting, numSort, new boolean[temp[0].getLength()], (int)(minPercent*sortTemp.length), sortByWeights,fullPFM,threshold,pivots);//TODO maxnum

		
		if(sortGlobally){
			
			ComparableElement<Integer, Double>[] pivAr = pivots.toArray( new ComparableElement[0] );
			
			Arrays.sort( pivAr );
			
			IntList sortPos = new IntList();
			DoubleList sortVals = new DoubleList();
			
			int off = 0;
			boolean[] exclude = new boolean[temp[0].getLength()];
			double[] vals = new double[exclude.length];
			for(int i=pivAr.length-1; i>=0 && off < numSort; i--){
				int idx = pivAr[i].getElement();
				double val = pivAr[i].getWeight();
				vals[idx] += val;
				/*if(!exclude[idx]){
					sortPos.add( idx );
					sortVals.add( val );
					exclude[idx] = true;
					off++;
				}*/
			}
			int[] o = ToolBox.order( vals, true );
			//o = new int[]{21,13,14};
			for(int i=0;i<o.length;i++){
				sortPos.add( o[i] );
				sortVals.add( vals[o[i]] );
			}
			
			sorted = partition( sortTemp, (int)(minPercent*sortTemp.length), sortPos.toArray(), sortVals.toArray(), 0, sortByWeights, null );
			
		}

		
		int totalHeight = chunkHeight;
		
		double[][] pwm = new double[seqs[0].getLength()][4];
		
		for(int i=0;i<sorted.length;i++){
			int numCurr = sorted[i].length;
			int heightCurr = totalHeight*numCurr/sortTemp.length;
			
			double meanW = 0.0;
			for(int k=0;k<numCurr;k++){
				meanW += sorted[i][k].getSecondElement();
			}
			meanW /= numCurr;
			for(int j=0;j<pwm.length;j++){
				Arrays.fill( pwm[j], 0.0 );
				for(int k=0;k<numCurr;k++){
					pwm[j][sorted[i][k].getFirstElement().discreteVal( j )]++;
				}
				Normalisation.sumNormalisation( pwm[j] );
			}
			
			plotColorLogo( graph, pwm, (meanW-mi)/(ma-mi), true, false, heightCurr, width, offx, offy );
			
			offy+=heightCurr;
			
		}
		
		return offx;
		
	}
	
	private static Pair<Sequence,Double>[] sortLocal3( Pair<Sequence, Double>[] sortTemp, int numBestForSorting, int maxNum, boolean[] exclude, int minElements, boolean sortByWeights, double[] prevInf ) {
		//exclude = exclude.clone();
		//System.out.println("maxNum: "+maxNum);
		
		if(maxNum == 0 || sortTemp.length < minElements*4){
			return sortTemp;
		}else{
			
			//double pc = 1.0;
			
			Sequence[] temp = new Sequence[sortTemp.length];
			for(int i=0;i<temp.length;i++){
				temp[i] = sortTemp[i].getFirstElement();
			}
			
			Pair<double[],double[][]> pair = getInformation( temp, numBestForSorting, exclude );
			double[] inf = pair.getFirstElement();
			int best = ToolBox.getMaxIndex( inf );
			
			
			
			int bestPrev = ToolBox.getMaxIndex( prevInf );
			
			//System.out.println(Arrays.toString( exclude ));
			//System.out.println(best+" "+inf[best]+" "+bestPrev+" "+prevInf[bestPrev]+" "+(temp.length+pc*16)+" "+numBestForSorting+" "+ToolBox.sum( exclude ));
			
			if(inf[best]/(temp.length)/(double)numBestForSorting < 0.1 && ( ToolBox.sum( exclude ) == 0 || prevInf[bestPrev]/(temp.length)/(double)ToolBox.sum( exclude ) < 0.1 )){
				return sortTemp;
			}else{
				
				if(inf[best]/(temp.length)/(double)numBestForSorting < 0.1){
					best = bestPrev;
				}
				
				double[] mis = pair.getSecondElement()[best];
				
				exclude[best] = true;
				//LinkedList<Sequence>[] partitions = new LinkedList[4];
				LinkedList<Pair<Sequence, Double>>[] partSort = new LinkedList[4];
				for(int i=0;i<partSort.length;i++){
					partSort[i] = new LinkedList<Pair<Sequence,Double>>();
				}
				
				double[] freq = new double[4];
				double[] ws = new double[4];
				for(int i=0;i<temp.length;i++){
					Sequence seq = temp[i];
					int idx = seq.discreteVal( best );
					//partitions[idx].add( seq );
					partSort[idx].add( sortTemp[i] );
					freq[idx]++;
					ws[idx] += sortTemp[i].getSecondElement(); 
				}
				if(sortByWeights){
					for(int i=0;i<freq.length;i++){
						freq[i] = ws[i]/freq[i];
					}
				}
				
				int[] ord = ToolBox.order( freq, true );
				//System.out.println(maxNum+": "+best);
				//System.out.println(Arrays.toString( freq ) );
				//System.out.println(Arrays.toString( ord ));
				int off = 0;
				for(int i=0;i<ord.length;i++){
					if(freq[ord[i]] > 0){
						
						double[] tempPrev = prevInf.clone();
						for(int j=0;j<tempPrev.length;j++){
							if(exclude[j]){
								tempPrev[j] = 0;
							}else{
								tempPrev[j] += mis[j];
								/*for(int k=0;k<depPos.length;k++){
									if(j==depPos[k]){
										tempPrev[j] = ToolBox.max( tempPrev[j], inf[j] );
									}
								}*/
							}
						}
						
						Pair<Sequence, Double>[] part = sortLocal3( partSort[ord[i]].toArray( new Pair[0] ), numBestForSorting, maxNum-1, exclude.clone(), minElements, sortByWeights, tempPrev );
						for(int j=0;j<part.length;j++,off++){
							sortTemp[off] = part[j];
						}
					}
				}
				return sortTemp;		
				
			}
			
			
		}
		
	}
	
	private static Pair<Sequence,Double>[][] sortLocal2( Pair<Sequence, Double>[] sortTemp, int numBestForSorting, int maxNum, boolean[] exclude, int minElements, boolean sortByWeights, double[][] fullPFM, double threshold, LinkedList<ComparableElement<Integer, Double>> pivots ) throws Exception {
		//exclude = exclude.clone();
		//System.out.println("maxNum: "+maxNum);
		
		if(maxNum == 0 || sortTemp.length < minElements){
			//System.out.println(maxNum+" "+sortTemp.length+" "+minElements);
			return new Pair[][]{sortTemp};
		}else{
						
			Sequence[] temp = new Sequence[sortTemp.length];
			for(int i=0;i<temp.length;i++){
				temp[i] = sortTemp[i].getFirstElement();
			}
			
			Pair<double[],double[][]> pair = getInformation( temp, numBestForSorting, exclude );
			double[] inf = pair.getFirstElement();
			double[][] deps = pair.getSecondElement();
			int best = ToolBox.getMaxIndex( inf );
			
			//System.out.println(Arrays.toString( inf ));
			//System.out.println(Arrays.toString( deps[best] ));
			//System.out.println(Arrays.toString( exclude ));
			//System.out.println(best+" "+inf[best]+" "+(temp.length)+" "+numBestForSorting);
			
			
			
			
			
			
			if(inf[best]/(double)(temp.length)/(double)numBestForSorting < threshold ){//TODO Threshold
				
				return new Pair[][]{sortTemp};
				
			}else{
								
				exclude[best] = true;
				
				
				Pair<Sequence, Double>[][] partSort = partition(sortTemp, minElements, best, inf[best]/(double)(temp.length)/(double)numBestForSorting, sortByWeights, pivots);

				
				double[] kls = deps[best];//getAvgKLs(sortTemp,partSort);
				for(int i=0;i<kls.length;i++){
					kls[i] /= (double)(temp.length);
				}
				//System.out.println("avg: "+Arrays.toString( kls ));
				
				int[] o2 = ToolBox.order( kls, true );
				int secpos = -1;
				for(int i=0;i<o2.length;i++){
					if(!exclude[o2[i]] && kls[o2[i]] > threshold){
						secpos = o2[i];
						break;
					}
				}
				if(secpos > -1){
					
					LinkedList<Pair<Sequence,Double>[]> list = new LinkedList<Pair<Sequence,Double>[]>();
					for(int i=0;i<partSort.length;i++){
						if(partSort[i] != null){
							Pair<Sequence, Double>[][] temp2 = partition( partSort[i], minElements, secpos, kls[secpos], sortByWeights, pivots );
							temp2 = joinSmall( temp2, minElements, sortByWeights );
							for(int j=0;j<temp2.length;j++){
								//if(temp2[j] != null){
									list.add( temp2[j] );
								//}
							}
						}
					}
					partSort = list.toArray( new Pair[0][0] );
					maxNum--;
					exclude[secpos] = true;
				}
				
				partSort = joinSmall(partSort,minElements,false);			
				
				
				LinkedList<Pair<Sequence,Double>[]> partitions = new LinkedList<Pair<Sequence,Double>[]>();
				for(int i=0;i<partSort.length;i++){
					if(partSort[i] != null){

						Pair<Sequence,Double>[][] part = sortLocal2( partSort[i], numBestForSorting, maxNum-1, exclude.clone(), minElements, sortByWeights,fullPFM, threshold, pivots );
						
						for(int j=0;j<part.length;j++){
							partitions.add( part[j] );
						}
					}
				}
				return partitions.toArray( new Pair[0][0] );	
				
			}
			
			
		}
		
	}
	
	
	private static Pair<Sequence, Double>[][] joinSmall( Pair<Sequence, Double>[][] partSort, int minElements, boolean sortByWeights ) {
		if(partSort.length == 1){
			return partSort;
		}
		
		boolean[] out = new boolean[partSort.length];
		int minAbove = Integer.MAX_VALUE;
		int idx = -1;
		int nout = 0;
		
		for(int i=0;i<partSort.length;i++){
			if(partSort[i] != null){
				if(partSort[i].length < minElements){
					out[i] = true;
					nout += partSort[i].length;
				}else{
					if(partSort[i].length < minAbove){
						minAbove = partSort[i].length;
						idx = i;
					}
				}
			}
		}
		
		if(nout == 0 && idx == -1){
			return partSort;
		}else if(idx == -1){
			minAbove = 0;
		}
		Pair<Sequence, Double>[] joined = new Pair[minAbove+nout];
		int off = 0;
		if(idx > -1){
			System.arraycopy( partSort[idx], 0, joined, 0, partSort[idx].length );
			off = partSort[idx].length;
		}else{
			for(int i=0;i<out.length;i++){
				if(out[i]){
					idx = i;
					break;
				}
			}
		}
		
		for(int i=0;i<partSort.length;i++){
			if(partSort[i] != null && out[i]){
				System.arraycopy( partSort[i], 0, joined, off, partSort[i].length );
				off += partSort[i].length;
				partSort[i] = null;
			}
		}
		partSort[idx] = joined;
		
		if(sortByWeights){
			double[] meanw = new double[partSort.length];
			for(int i=0;i<partSort.length;i++){
				if(partSort[i] == null){
					meanw[i] = Double.NEGATIVE_INFINITY;
				}else{
					for(int j=0;j<partSort[i].length;j++){
						meanw[i] += partSort[i][j].getSecondElement();
					}
					meanw[i] /= partSort[i].length;
				}
			}
				
			int[] o = ToolBox.order( meanw, true );
			
			Pair<Sequence, Double>[][] temp = new Pair[partSort.length][];
			for(int i=0;i<o.length;i++){
				temp[i] = partSort[o[i]];
			}
			
			partSort = temp;
			
		}else{//TODO FIXME remove comment ???

			/*double[] nums = new double[partSort.length];
			for(int i=0;i<nums.length;i++){
				if(partSort[i] != null){
					nums[i] = partSort[i].length;
				}else{
					nums[i] = 0;
				}
			}
			int[] o = ToolBox.order( nums, true );
			
			Pair<Sequence, Double>[][] temp = new Pair[partSort.length][];
			for(int i=0;i<o.length;i++){
				temp[i] = partSort[o[i]];
			}
			
			partSort = temp;
			*/
		}
		
		return partSort;
		
	}

	private static double[] getAvgKLs( Pair<Sequence, Double>[] sortTemp, Pair<Sequence, Double>[][] partSort ) throws EmptyDataSetException, WrongAlphabetException, CloneNotSupportedException {
		
		double[][] all = getPFM(sortTemp);
		double n = sortTemp.length;
		
		double[] allAvgKLs = new double[sortTemp[0].getFirstElement().getLength()];
		
		for(int i=0;i<partSort.length;i++){
			if(partSort[i] != null){
				double[][] curr = getPFM(partSort[i]);
				double m = partSort[i].length;
				double[] kls = getKLDivergence( curr, ArrayHandler.clone( all ) );
				for(int j=0;j<kls.length;j++){
					allAvgKLs[j] += m/n* kls[j];
				}
			}
		}
		return allAvgKLs;
		
	}

	private static double[][] getPFM( Pair<Sequence, Double>[] sortTemp ) throws EmptyDataSetException, WrongAlphabetException {
		LinkedList<Sequence> seqTemp = new LinkedList<Sequence>();
		for(int i=0;i<sortTemp.length;i++){
			seqTemp.add( sortTemp[i].getFirstElement() );
		}
		
		
		double[][] partPFM = PFMComparator.getPFM( new DataSet("",seqTemp) );
		return partPFM;
	}

	private static Pair<Sequence, Double>[][] partition( Pair<Sequence,Double>[] sortTemp, int minElements, int[] ord, double[] values, int idx, boolean sortByWeights, LinkedList<ComparableElement<Integer, Double>> pivots ){
		
		if(idx == ord.length){
			return new Pair[][]{sortTemp};
		}else{
			
			
			Pair<Sequence, Double>[][] partitions = partition(sortTemp,minElements,ord[idx],values[idx],sortByWeights,pivots);
			
			LinkedList<Pair<Sequence,Double>[]> pairs = new LinkedList<Pair<Sequence,Double>[]>();
			for(int i=0;i<partitions.length;i++){
				if(partitions[i] != null){
					Pair<Sequence,Double>[][] part2 = partition( partitions[i],minElements,ord,values,idx+1,sortByWeights, pivots );
					for(int j=0;j<part2.length;j++){
						pairs.add( part2[j] );
					}
				}
			}
			
			return pairs.toArray( new Pair[0][0] );
			
		}
		
	}

	private static Pair<Sequence, Double>[][] partition( Pair<Sequence, Double>[] sortTemp, int minElements,int curr, double value, boolean sortByWeights, LinkedList<ComparableElement<Integer, Double>> pivots ) {
			
		if(sortTemp.length < minElements){
			return new Pair[][]{sortTemp};
		}
		
			if(pivots != null){
				pivots.add( new ComparableElement<Integer, Double>( curr, value*sortTemp.length ) );
			}
		
			LinkedList<Pair<Sequence, Double>>[] partSort = new LinkedList[4];
			for(int i=0;i<partSort.length;i++){
				partSort[i] = new LinkedList<Pair<Sequence,Double>>();
			}
			
			double[] freq = new double[4];
			double[] ws = new double[4];
			for(int i=0;i<sortTemp.length;i++){
				Sequence seq = sortTemp[i].getFirstElement();
				int idx = seq.discreteVal( curr );
				//partitions[idx].add( seq );
				partSort[idx].add( sortTemp[i] );
				freq[idx]++;
				ws[idx] += sortTemp[i].getSecondElement(); 
			}
			if(sortByWeights){
				for(int i=0;i<freq.length;i++){
					freq[i] = ws[i]/freq[i];
				}
			}
			
			int[] ord = ToolBox.order( freq, true );
			//System.out.println(maxNum+": "+best);
			//System.out.println(Arrays.toString( freq ) );
			//System.out.println(Arrays.toString( ord ));
			int off = 0;
			Pair<Sequence,Double>[][] partitions = new Pair[4][];
			for(int i=0;i<ord.length;i++){
				if(partSort[ord[i]].size() > 0){
					
					partitions[i] = partSort[ord[i]].toArray( new Pair[0] );
							//sortLocal2( , numBestForSorting, maxNum-1, exclude.clone(), minElements, sortByWeights,fullPFM );
				}
			}
			return partitions;
		
	}

	private static double[] getKLDivergence( double[][] partPFM, double[][] fullPFM ) {
		
		double full = ToolBox.sum( fullPFM[0] );
		double part = ToolBox.sum( partPFM[0] );
		
		
		double[] kls = new double[partPFM.length];
		
		for(int i=0;i<fullPFM.length;i++){
			Normalisation.sumNormalisation( fullPFM[i] );
			Normalisation.sumNormalisation( partPFM[i] );
			
			for(int j=0;j<fullPFM[i].length;j++){
				if(partPFM[i][j] > 0){
					kls[i] += partPFM[i][j] * Math.log( partPFM[i][j] / fullPFM[i][j] );
				}
				
			}
			
			
		}
		
		return kls;
		
	}

	/*public static int[] getOrder( Sequence[] temp, int numBestForSorting, int numSortingPositions ) {
		
		IntList sortPos = new IntList();
		boolean[] exclude = new boolean[temp[0].getLength()];
		
		
		double[] inf = getInformation( temp, numBestForSorting, exclude ).getFirstElement();
		int best = ToolBox.getMaxIndex( inf );
		exclude[best] = true;
		if(inf[best]/(temp.length)/(double)numBestForSorting > 0.1){
			sortPos.add( best );
		}
		//System.out.println("best 0: "+best+" "+inf[best]/(temp.length+pc*16)/numBestForSorting);
		for(int i=1;i<numSortingPositions;i++){
			//pc /= 4.0;
			
			inf = getInformation( temp, numBestForSorting, exclude ).getFirstElement();
			best = ToolBox.getMaxIndex( inf );
			exclude[best] = true;
			//System.out.println("best "+i+": "+best+" "+inf[best]/(double)(temp.length+pc*16)/(double)numBestForSorting);
			if(inf[best]/(double)(temp.length)/(double)numBestForSorting > 0.1){
				sortPos.add( best );
			}
		}
		return sortPos.toArray();
	}*/
	
	public static LinkedList<Sequence>[] getPartitions(Sequence[] temp, IntList sortPos){
		LinkedList<Sequence>[] lists = new LinkedList[(int)Math.pow( 4, sortPos.length() )];
		for(int i=0;i<lists.length;i++){
			lists[i] = new LinkedList<Sequence>();
		}
		int[] pows = new int[sortPos.length()];
		for(int i=0;i<pows.length;i++){
			pows[i] = (int) Math.pow( 4, i );
		}
		for(int i=0;i<temp.length;i++){
			int idx = 0;
			for(int j=0;j<sortPos.length();j++){
				idx += pows[j]*temp[i].discreteVal( sortPos.get( j ) );
			}
			lists[idx].add( temp[i] );
		}
		return lists;
	}

	public static Pair<double[],double[][]> getInformation(Sequence[] seqs, int numBest, boolean[] exclude){
			
		
		double[][] mis = getDeps( seqs, null );
		double[] vals = new double[seqs[0].getLength()];
		
		//int[][] depPoss = new int[vals.length][0];
		for(int i=0;i<vals.length;i++){
			//Arrays.sort( mis[i] );
			int[] o = ToolBox.order( mis[i], true );
			vals[i] = 0;
			//IntList depPos = new IntList();
			if(exclude == null || /*!exclude[o[i]]*/!exclude[i]){
				for(int j=0,k=0;k<numBest && j<mis[i].length;j++){
					//vals[i] += mis[i][mis[i].length-j-1];
					//if(exclude==null || !exclude[o[j]]){//TODO start exclude
						vals[i] += mis[i][o[j]];
						k++;
					//}//TODO end exclude
					//depPos.add( o[j] );
				}
			}
			//depPoss[i] = depPos.toArray();
		}
		//System.out.println(numBest+" "+Arrays.toString( vals ));
		return new Pair<double[],double[][]>(vals,mis);
	}
	
	public static void plotColorLogo(Graphics2D graph, double[][] pwm, double weight, boolean mix, boolean icscale, int height, int width, int offx, int offy){
		graph = (Graphics2D)graph.create();
		
		Color back = graph.getColor();
		graph.setColor( Color.WHITE );
		graph.fillRect( offx, offy, width, height );
		graph.setColor( back );
		
		int partWidth = (int)Math.floor( (double)width/(double)(pwm.length) );
				
		double[] temp = new double[4];
		
		for(int i=0;i<pwm.length;i++){
			
			if(mix){
				Color c = getColor( pwm[i], icscale ? (float)Math.sqrt( getICScale( pwm[i] ) ) : 1f );
				back = graph.getColor();
				graph.setColor( c );
				graph.fillRect( offx, offy, partWidth, height );
				graph.setColor( back );
			}else{
				int off2 = 0;
				float ic = icscale ? (float)Math.sqrt( getICScale( pwm[i] ) ) : 1f;
				//System.out.println(ic);
				for(int j=0;j<pwm[i].length;j++){
					Arrays.fill( temp, 0 );
					temp[j] = 1.0;
					Color c = getColor( temp, ic );
					back = graph.getColor();
					graph.setColor( c );
					int partHeight = (int)Math.round( height*pwm[i][j] );
					graph.fillRect( offx, offy+off2, partWidth, j<pwm[i].length-1 ? partHeight : height-off2 );
					off2 += partHeight;
					graph.setColor( back );
				}
			}
			
			offx += partWidth;
			
		}
		
		graph.setColor( back );
		
	}
	
	/**
	 * Plots the sequence logo for the position weight matrix given in <code>ps</code>. 
	 * The rows of <code>ps</code> correspond to positions in the sequence logo. Each row must be normalized.
	 * For a given <code>height</code> (in pixels), the width is chosen automatically depending on the number of rows
	 * in <code>ps</code>.
	 * 
	 * The positions of the sequence logo are numbered contiguously from 1 to <code>ps.length</code>. The label of the
	 * x-axis is set to &quot;Position&quot;, and the label of the y-axis is set to &quot;bits&quot;.
	 * 
	 * The sequence logo is written to the PNG file given in <code>path</code>.
	 * 
	 * @param path the path to the PNG file written
	 * @param height the height of the PNG image (in pixels)
	 * @param ps the position weight matrix
	 * @throws IOException if the file could not be written
	 */
	public static void plotLogoToPNG(String path, int height, double[][] ps) throws IOException{
		Pair<BufferedImage, Graphics2D> pair = getBufferedImageAndGraphics( height, ps );
		Graphics2D g = pair.getSecondElement();
		g.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
		g.setRenderingHint(RenderingHints.KEY_RENDERING, RenderingHints.VALUE_RENDER_QUALITY);
		plotLogo( g, height, ps );
		ImageIO.write( pair.getFirstElement(), "png", new File(path) );
	}
	
	/**
	 * Plots the sequence logo for the position weight matrix given in <code>ps</code>. 
	 * The rows of <code>ps</code> correspond to positions in the sequence logo. Each row must be normalized.
	 * For a given <code>height</code> (in pixels), the width is chosen automatically depending on the number of rows
	 * in <code>ps</code>.
	 * 
	 * The positions of the sequence logo are numbered contiguously from 1 to <code>ps.length</code>. The label of the
	 * x-axis is set to &quot;Position&quot;, and the label of the y-axis is set to &quot;bits&quot;.
	 * 
	 * The sequence logo is return as a {@link BufferedImage}.
	 * 
	 * @param height the height of the PNG image (in pixels)
	 * @param ps the position weight matrix
	 * @return the sequence logo
	 */
	public static BufferedImage plotLogoToBufferedImage(int height, double[][] ps) {
		Pair<BufferedImage, Graphics2D> pair = getBufferedImageAndGraphics( height, ps );
		Graphics2D g = pair.getSecondElement();
		g.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
		g.setRenderingHint(RenderingHints.KEY_RENDERING, RenderingHints.VALUE_RENDER_QUALITY);
		plotLogo( g, height, ps );
		return pair.getFirstElement();
	}
	
	/**
	 * Creates a new {@link BufferedImage} with given height and width chosen automatically according to the number of rows
	 * of <code>ps</code>, and returns this {@link BufferedImage} and its {@link Graphics2D} object.
	 * 
	 * @param height the height (in pixels)
	 * @param ps the position weight matrix
	 * @return the created {@link BufferedImage} and its {@link Graphics2D} object
	 */
	public static Pair<BufferedImage,Graphics2D> getBufferedImageAndGraphics(int height,double[][] ps){
		int w = getWidth( height, ps );
		BufferedImage img = new BufferedImage( w, height, BufferedImage.TYPE_INT_RGB);
		Graphics2D g = (Graphics2D)img.getGraphics();
		
		return new Pair<BufferedImage, Graphics2D>( img, g );
	}
	
	/**
	 * Plots the sequence logo for the position weight matrix given in <code>ps</code>. 
	 * The rows of <code>ps</code> correspond to positions in the sequence logo. Each row must be normalized.
	 * For a given <code>height</code> (in pixels), the width is chosen automatically depending on the number of rows
	 * in <code>ps</code>.
	 * 
	 * The positions of the sequence logo are numbered contiguously from 1 to <code>ps.length</code>. The label of the
	 * x-axis is set to &quot;Position&quot;, and the label of the y-axis is set to &quot;bits&quot;.
	 * 
	 * The sequence logo is written to the {@link Graphics2D} object given in <code>g</code>.
	 * 
	 * @param g the {@link Graphics2D} object
	 * @param h the height of the sequence logo
	 * @param ps the position weight matrix
	 */
	public static void plotLogo(Graphics2D g, int h, double[][] ps){
		plotLogo( g, h, ps, null, "Position", "bits" );
	}
	
	/**
	 * Plots the sequence logo for the position weight matrix given in <code>ps</code>. 
	 * The rows of <code>ps</code> correspond to positions in the sequence logo. Each row must be normalized.
	 * For a given <code>height</code> (in pixels), the width is chosen automatically depending on the number of rows
	 * in <code>ps</code>.
	 * 
	 * The sequence logo is written to the {@link Graphics2D} object given in <code>g</code>.
	 * 
	 * @param g the {@link Graphics2D} object
	 * @param height the height of the sequence logo
	 * @param ps the position weight matrix
	 * @param labels the labels of the positions of the sequence logo, if <code>null</code> the positions are numbered contiguously from 1 to <code>ps.length</code>
	 * @param labX the label of the x-axis
	 * @param labY the label of the y-axis
	 */
	public static void plotLogo(Graphics2D g, int height, double[][] ps, String[] labels, String labX, String labY){
		int w = getWidth( height, ps );
		plotLogo( g, w, height, ps, labels, labX, labY );
	}
	
	/**
	 * Returns the automatically chosen width for a given height and position weight matrix.
	 * @param height the height
	 * @param ps the position weight matrix
	 * @return the width
	 */
	public static int getWidth(int height, double[][] ps){
		return getWidth(height, ps.length);
	}
	
	public static int getWidth(int height, int numCol){
		return (int)(height/6.0*(numCol+1.5));
	}
	
	public static int getColumnWidth(int height){
		return (int)(height/6.0);
	}
	
	/**
	 * Returns the automatically chosen height for a given width and position weight matrix.
	 * @param width the width
	 * @param ps the position weight matrix
	 * @return the height
	 */
	public static int getHeight(int width, double[][] ps){
		return (int)(width*6.0/(ps.length+1.5));
	}
	
	/**
	 * Plots the sequence logo for the position weight matrix given in <code>ps</code>. 
	 * The rows of <code>ps</code> correspond to positions in the sequence logo. Each row must be normalized.
	 * 
	 * The sequence logo is written to the {@link Graphics2D} object given in <code>g</code>.
	 * 
	 * @param g the {@link Graphics2D} object
	 * @param w the width of the sequence logo
	 * @param h the height of the sequence logo
	 * @param ps the position weight matrix
	 * @param labels the labels of the positions of the sequence logo, if <code>null</code> the positions are numbered contiguously from 1 to <code>ps.length</code>
	 * @param labX the label of the x-axis
	 * @param labY the label of the y-axis
	 */
	public static void plotLogo(Graphics2D g, int w, int h, double[][] ps, String[] labels, String labX, String labY){
		plotLogo( g, 0, h, w, h, ps, labels, labX, labY );
	}
	
	/**
	 * Plots the sequence logo for the position weight matrix given in <code>ps</code>. 
	 * The rows of <code>ps</code> correspond to positions in the sequence logo. Each row must be normalized.
	 * 
	 * The sequence logo is written to the {@link Graphics2D} object given in <code>g</code>.
	 * 
	 * @param g the {@link Graphics2D} object
	 * @param x the x-coordinate of the bottom left corner of the sequence logo
	 * @param y the y-coordinate of the bottom left corner of the sequence logo (<code>-h</code> results in a sequence logo spanning from <code>0</code> to <code>h</code>)
	 * @param w the width of the sequence logo
	 * @param h the height of the sequence logo
	 * @param ps the position weight matrix
	 * @param labels the labels of the positions of the sequence logo, if <code>null</code> the positions are numbered contiguously from 1 to <code>ps.length</code>
	 * @param labX the label of the x-axis
	 * @param labY the label of the y-axis
	 */
	public static void plotLogo(Graphics2D g, int x, int y, int w, int h, double[][] ps, String[] labels, String labX, String labY){
		g = (Graphics2D)g.create();
		/*g.scale( 1.0/(h*100), 1.0/(h*100) );
		x *= h*100;
		y *= h*100;
		w *= h*100;
		h *= h*100;*/
		g.setColor( Color.WHITE );
		g.fillRect( x, y-h, w, h );
		
		Font font = new Font(g.getFont().getName(),Font.BOLD,h/10);//17
		g.setFont( font );
		
		if(labels == null){
			labels = new String[ps.length];
			for(int i=0;i<ps.length;i++){
				labels[i] = (i+1)+"";
			}
		}
		

		double wl = h*0.4;
		
		double w2 = (w-wl)/(ps.length);
		double x2 = wl*0.9;
		
		double h2 = h*0.65;
		double y2 = y - h*0.3;
		
		g.setColor( Color.BLACK );
		g.setStroke( new BasicStroke( h/100f ) );
		g.drawLine( x+(int)x2, (int)(y2+0.05*h), x+(int)(x2+w2*ps.length), (int)(y2+0.05*h) );
		g.drawLine( x+(int)(x2*0.94), (int)y2, x+(int)(x2*0.94), (int)(y2-h2) );
		String[] labs = {"0", "0.5", "1", "1.5", "2"};
		for(int i=0;i<=4;i++){
			g.drawLine( x+(int)(x2*0.7), (int)(y2-i*h2/4.0), x+(int)(x2*0.94), (int)(y2-i*h2/4.0) );
			Rectangle2D rect = g.getFontMetrics().getStringBounds( labs[i], g );
			g.drawString( labs[i], x+(int)(x2*0.6-rect.getWidth()), (int)(y2-i*h2/4.0 - rect.getCenterY()) );
		}
		AffineTransform back = g.getTransform();
		g.rotate( -Math.PI/2 );
		Rectangle2D rect = g.getFontMetrics().getStringBounds( labY, g );

		g.drawString(labY,-(int)(y2-2*h2/4.0 + rect.getCenterX()), (int)(x+rect.getWidth()/2d));
		g.setTransform( back );
		
		rect = g.getFontMetrics().getStringBounds( labX, g );
		g.drawString( labX, x+(int)(x2+w2*ps.length/2.0-rect.getCenterX()), (int)(y-0.2*rect.getHeight()) );
		for(int i=0;i<ps.length;i++){
			plotLogo( g, x+x2, y2, w2, h2, ps[i] );
			g.setColor( Color.BLACK );
			rect = g.getFontMetrics().getStringBounds( labels[i], g );
			//g.drawLine( (int)(x2+w2/2d), (int)(y2+0.05*h), (int)(x2+w2/2d), (int)(y2+0.05*h+0.05*h ) );
			g.drawString( labels[i], x+(float)(x2+w2/2d-rect.getCenterX()), (float)(y2+1.5*rect.getHeight()) );
			x2 += w2;
		}
	}
	
	/**
	 * Plots the TALgetter logo for the binding specificities given in <code>ps</code>. 
	 * The rows of <code>ps</code> correspond to specificities of the RVDs. Each row must be normalized.
	 * For a given <code>height</code> (in pixels), the width is chosen automatically depending on the number of rows
	 * in <code>ps</code>. In addition, the importance of RVDs is plotted as given in <code>imp</code>.
	 * 
	 * The labels of the RVDs are given in <code>lab</code>. The label of the
	 * x-axis is set to &quot;RVD&quot;, and the label of the y-axes are set to &quot;bits&quot; and &quot;Importance&quot;, respectively.
	 * 
	 * The TALgetter logo is written to the PNG file given in <code>path</code>.
	 * 
	 * @param path the path to the PNG file written
	 * @param height the height of the PNG image (in pixels)
	 * @param ps the binding specificities of RVDs
	 * @param imp the importance of RVDs
	 * @param lab the amino acids of the RVDs in one-letter code
	 * @throws IOException if the file could not be written
	 */
	public static void plotTALgetterLogoToPNG(String path, int height, double[][] ps, double[] imp, String[] lab) throws IOException{
		Pair<BufferedImage, Graphics2D> pair = getBufferedImageAndGraphics( height, ps );
		Graphics2D g = pair.getSecondElement();
		g.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
		g.setRenderingHint(RenderingHints.KEY_RENDERING, RenderingHints.VALUE_RENDER_QUALITY);
		int w = getWidth( height, ps );
		
		plotTALgetterLogo( g, 0, height, w, height, ps, imp, lab, "RVD", "bits", "Importance" );
		ImageIO.write( pair.getFirstElement(), "png", new File(path) );
	}
	
	/**
	 * Plots the TALgetter logo for the binding specificities given in <code>ps</code>. 
	 * The rows of <code>ps</code> correspond to specificities of the RVDs. Each row must be normalized.
	 * For a given <code>height</code> (in pixels), the width is chosen automatically depending on the number of rows
	 * in <code>ps</code>. In addition, the importance of RVDs is plotted as given in <code>imp</code>.
	 * 
	 * The labels of the RVDs are given in <code>lab</code>. The label of the
	 * x-axis is set to &quot;RVD&quot;, and the label of the y-axes are set to &quot;bits&quot; and &quot;Importance&quot;, respectively.
	 * 
	 * The TALgetter logo is returned as {@link BufferedImage}.
	 * 
	 * @param path the path to the PNG file written
	 * @param height the height of the PNG image (in pixels)
	 * @param ps the binding specificities of RVDs
	 * @param imp the importance of RVDs
	 * @param lab the amino acids of the RVDs in one-letter code
	 * @return the TALgetter logo
	 */
	public static BufferedImage plotTALgetterLogoToBufferedImage(int height, double[][] ps, double[] imp, String[] lab) {
		Pair<BufferedImage, Graphics2D> pair = getBufferedImageAndGraphics( height, ps );
		Graphics2D g = pair.getSecondElement();
		g.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
		g.setRenderingHint(RenderingHints.KEY_RENDERING, RenderingHints.VALUE_RENDER_QUALITY);
		int w = getWidth( height, ps );
		
		plotTALgetterLogo( g, 0, height, w, height, ps, imp, lab, "RVD", "bits", "Importance" );
		return pair.getFirstElement();
	}
	
	/**
	 * Plots the TALgetter logo for the binding specificities given in <code>ps</code>. 
	 * The rows of <code>ps</code> correspond to specificities of the RVDs. Each row must be normalized.
	 * For a given <code>height</code> (in pixels), the width is chosen automatically depending on the number of rows
	 * in <code>ps</code>. In addition, the importance of RVDs is plotted as given in <code>imp</code>.
	 * 
	 * The labels of the RVDs are given in <code>lab</code>. The label of the
	 * x-axis is set to &quot;RVD&quot;, and the label of the y-axes are set to &quot;bits&quot; and &quot;Importance&quot;, respectively.
	 * 
	 * The TALgetter logo is returned as {@link BufferedImage}.
	 * 
	 * @param g the {@link Graphics2D} object
	 * @param x the x-coordinate of the bottom left corner of the TALgetter logo
	 * @param y the y-coordinate of the bottom left corner of the TALgetter logo (<code>-h</code> results in a logo spanning from <code>0</code> to <code>h</code>)
	 * @param w the width of the TALgetter logo
	 * @param h the height of the TALgetter logo
	 * @param ps the binding specificities of RVDs
	 * @param imp the importance of RVDs
	 * @param labels the amino acids of the RVDs in one-letter code
	 * @param labX the label of the x-axis
	 * @param labY the label of the y-axis
	 * @param labY2 the label of the second (importance) y-axis
	 */
	public static void plotTALgetterLogo(Graphics2D g, int x, int y, int w, int h, double[][] ps, double[] imp, String[] labels, String labX, String labY, String labY2){
		g = (Graphics2D)g.create();
		g.setColor( Color.WHITE );
		g.fillRect( x, y-h, w, h );
		

		Font font = new Font(g.getFont().getName(),Font.BOLD,h/17);
		g.setFont( font );
		
		if(labels == null){
			labels = new String[ps.length];
			for(int i=0;i<ps.length;i++){
				labels[i] = (i+1)+"";
			}
		}
		

		double wl = h*0.5;
		
		double w2 = (w-wl)/(ps.length);
		double x2 = x + wl*0.45;
		
		double h2 = h*0.7;
		double y2 = y*0.75;
		
		g.setColor( Color.BLACK );
		g.setStroke( new BasicStroke( h/400+1 ) );
		g.drawLine( (int)x2, (int)(y2*1.04)+1, (int)(x2+w2*ps.length), (int)(y2*1.04)+1 );
		g.drawLine( (int)(x2*0.94)-1, (int)y2, (int)(x2*0.94)-1, (int)(y2-h2) );
		String[] labs = {"0", "0.5", "1", "1.5", "2"};
		for(int i=0;i<=4;i++){
			g.drawLine( (int)(x2*0.7), (int)(y2-i*h2/4.0), (int)(x2*0.94)-1, (int)(y2-i*h2/4.0) );
			Rectangle2D rect = g.getFontMetrics().getStringBounds( labs[i], g );
			g.drawString( labs[i], (int)(x2*0.6-rect.getWidth()-2), (int)(y2-i*h2/4.0 - rect.getCenterY()) );
		}
		AffineTransform back = g.getTransform();
		g.rotate( -Math.PI/2 );
		Rectangle2D rect = g.getFontMetrics().getStringBounds( labY, g );

		g.drawString(labY,-(int)(y2-2*h2/4.0 + rect.getCenterX()), (int)(x+rect.getHeight()));
		g.setTransform( back );
		
		
		rect = g.getFontMetrics().getStringBounds( labX, g );
		g.drawString( labX, (int)(x2+w2*ps.length/2.0-rect.getCenterX()), (int)(y-0.3*rect.getHeight()) );
		for(int i=0;i<ps.length;i++){
			plotLogo( g, x2, y2, w2, h2, ps[i] );
			g.setColor( Color.BLACK );
			rect = g.getFontMetrics().getStringBounds( labels[i], g );
			g.drawString( labels[i], (float)(x2+w2/2d-rect.getCenterX()), (float)(y2+2*rect.getHeight()) );
			x2 += w2;
		}
		
		g.drawLine( (int)(x2+w2*0.1)+1, (int)y2, (int)(x2+w2*0.1)+1, (int)(y2-h2) );
		labs = new String[]{"0", "0.5", "1"};
		Rectangle2D rect2 = g.getFontMetrics().getStringBounds( "0.5", g );
		for(int i=0;i<=2;i++){
			g.drawLine( (int)(x2+w2*0.34)+2, (int)(y2-i*h2/2.0), (int)(x2+w2*0.1)+1, (int)(y2-i*h2/2.0) );
			//rect = g.getFontMetrics().getStringBounds( labs[i], g );
			g.drawString( labs[i], (int)((x2+w2*0.5)+2), (int)(y2-i*h2/2.0 - rect.getCenterY()) );
		}
		
		back = g.getTransform();
		g.rotate( -Math.PI/2 );
		rect = g.getFontMetrics().getStringBounds( labY2, g );

		g.drawString(labY2,-(int)(y2-2*h2/4.0 + rect.getCenterX()), (int)(w-rect.getHeight()/2d));
		g.setTransform( back );
		
		x2 = x + wl*0.45 + 2*w2;
		

		g.setColor( Color.GRAY );
		for(int i=1;i<imp.length;i++){
			g.drawLine( (int)(x2-w2/2d + w2/20d),(int)(y2-h2*imp[i-1]+w2/20d) , (int)(x2+w2/2d - w2/20d), (int)(y2-h2*imp[i]+w2/20d) );
			x2 += w2;
		}
		
		x2 = x + wl*0.45 + w2;
		
		g.setColor( Color.BLUE );
		for(int i=0;i<imp.length;i++){
			
			g.fillRect( (int)(x2+w2/2d - w2/20d), (int)(y2-h2*imp[i]) , (int)(w2/10d), (int)(w2/10d) );
			x2 += w2;
		}
		
	}
	
	protected static void plotLogo(Graphics2D g, double x, double y, double w, double h, double[] p){
		
		//y += h;
		
		double ic = getICScale(p);
		h *= ic;
		//System.out.println("h: "+h);
		
		double[] mp = p.clone();
		for(int i=0;i<mp.length;i++){
			mp[i] *= -1;
		}
		
		int[] r = ToolBox.rank( mp, false );
		int[] order = new int[r.length];
		for(int i=0;i<r.length;i++){
			order[r[i]] = i;
		}
		
		
		for(int i=0;i<order.length;i++){
			double curr = p[order[i]];
			if(order[i] == 0){
				g.setColor( Color.GREEN );
				g.fill(getA( x, y, w, h*curr ));
			}else if(order[i] == 1){
				g.setColor( Color.BLUE );
				g.fill(getC( x, y, w, h*curr ));
			}else if(order[i] == 2){
				g.setColor( Color.ORANGE );
				g.fill(getG( x, y, w, h*curr ));
			}else{
				g.setColor( Color.RED );
				g.fill(getT( x, y, w, h*curr ));
			}
			//System.out.println("y: "+y);
			y -= h*curr;
		}
		//System.out.println(y);
	}
	
	public static double getICScale( double[] p ) {
		double ic = Math.log( 4 )/Math.log( 2 );
		for(int i=0;i<p.length;i++){
			if(p[i] > 0){
				ic += p[i]*Math.log( p[i] )/Math.log( 2 );
			}
		}
		ic /= 2.0;
		return ic;
	}

	private static Area getC( double x, double y, double w, double h ){
		
		Shape s = new Ellipse2D.Double( 0, -90, 90, 90 );
		Area a1 = new Area( s );
		
		Shape s2 = new Ellipse2D.Double( 15, -75, 60, 60 );
		Area a2 = new Area(s2);
		
		a1.subtract( a2 );
		
		Shape s3 = new Rectangle2D.Double( 65, -60, 30, 30 );
		Area a3 = new Area( s3 );
		
		a1.subtract( a3 );
		
		AffineTransform t = new AffineTransform();
		t.scale( 1d/88d, 1d/90d );
		t.scale( w, h );
		a1.transform( t );
		t = new AffineTransform();
		t.translate( x, y );
		a1.transform( t );
		return a1;
	}
	
	private static Area getT( double x, double y, double w, double h ){
		
		Shape s = new Rectangle2D.Double( 37.5, -100, 15, 100 );
		Area a1 = new Area( s );
		
		Shape s2 = new Rectangle2D.Double(0,-100,90,15);
		Area a2 = new Area( s2 );
		
		a1.add( a2 );
		
		AffineTransform t = new AffineTransform();
		t.scale( 1d/90d, 1d/100d );
		t.scale( w, h );
		a1.transform( t );
		t = new AffineTransform();
		t.translate( x, y );
		a1.transform( t );
		return a1;
		
	}
	
	private static Area getG( double x, double y, double w, double h ){
		
		Shape s = new Ellipse2D.Double( 0, -90, 90, 90 );
		Area a1 = new Area( s );
		
		Shape s2 = new Ellipse2D.Double( 15, -75, 60, 60 );
		Area a2 = new Area(s2);
		
		a1.subtract( a2 );
		
		Shape s3 = new Rectangle2D.Double( 65, -60, 30, 30 );
		Area a3 = new Area( s3 );
		
		a1.subtract( a3 );
		
		Shape s4 = new Rectangle2D.Double( 55, -40, 35, 15 );
		Area a4 = new Area( s4 );
		
		a1.add( a4 );
		
		Shape s5 = new Rectangle2D.Double( 80, -40, 10, 40 );
		Area a5 = new Area( s5 );
		
		a1.add( a5 );
		
		AffineTransform t = new AffineTransform();
		t.scale( 1d/90d, 1d/90d );
		t.scale( w, h );
		a1.transform( t );
		t = new AffineTransform();
		t.translate( x, y );
		a1.transform( t );
		return a1;
	}
	
	private static Area getA( double x, double y, double w, double h ){
		
		Shape s = new Polygon( new int[]{0,40,50,90,75,45,45,15,0}, new int[]{0,-100,-100,0,0,-80,-80,0,0}, 9 );
		Area a = new Area( s );
		
		Shape s2 = new Polygon(new int[]{20,70,70,20},new int[]{-35,-35,-50,-50},4);
		Area a2 = new Area( s2 );
		
		a.add( a2 );
		
		AffineTransform t = new AffineTransform();
		t.scale( 1.0/90.0, 1.0/100.0 );
		t.scale( w, h );
		a.transform( t );
		t = new AffineTransform();
		t.translate( x, y );
		a.transform( t );
		return a;
	}
	
	
	
}
