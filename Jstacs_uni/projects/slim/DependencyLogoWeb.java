package projects.slim;
import java.awt.Graphics2D;
import java.awt.RenderingHints;
import java.awt.image.BufferedImage;

import projects.dimont.DimontWeb;
import de.jstacs.data.DataSet;
import de.jstacs.io.FileManager;
import de.jstacs.utils.Pair;
import de.jstacs.utils.SeqLogoPlotter;
import de.jstacs.utils.galaxy.GalaxyAdaptor;
import de.jstacs.utils.galaxy.GalaxyAdaptor.FileResult;
import de.jstacs.utils.galaxy.GalaxyAdaptor.LinkedImageResult;
import de.jstacs.utils.galaxy.GalaxyAdaptor.Protocol;
import de.jstacs.utils.graphics.GraphicsAdaptor;
import de.jstacs.utils.graphics.GraphicsAdaptorFactory;
import de.jstacs.utils.graphics.GraphicsAdaptorFactory.OutputFormat;
import de.jstacs.utils.graphics.RasterizedAdaptor;



public class DependencyLogoWeb {

	
	/**
	 * @param args
	 */
	public static void main( String[] args ) throws Exception {
		
		DependencyLogoWebParameterSet params = new DependencyLogoWebParameterSet();
		
		boolean[] lines = new boolean[params.getNumberOfParameters()];
		
		GalaxyAdaptor ga = new GalaxyAdaptor( params,lines,"Dependency logo", " Plot dependency logos from a tabular file", "0.1", "java -Xms256M -Xmx2G -Djava.awt.headless=true -jar "+System.getProperty( "user.dir" )+System.getProperty( "file.separator" )+"DependencyLogoWeb.jar", "jobname" );
		ga.setHelp( FileManager.readInputStream( DimontWeb.class.getClassLoader().getResourceAsStream( "projects/slim/helpLogo.txt" ) ).toString() );//TODO
		
		if(!ga.parse( args, false )){
			System.exit( 1 );
		}	
		
		Protocol prot = ga.getProtocol( false );
		
		LinkedImageResult res = run(params, "./dependency_logo", prot);
		
		if(res != null){
			ga.addResult( res, true, true, res.getLink().getExtension() );
		}

		ga.writeOutput();
		
	}
	
	public static LinkedImageResult run(DependencyLogoWebParameterSet params, String pdfPath, Protocol prot) throws Exception{
		
		Pair<DataSet,double[]> pair = params.getData();
		
		OutputFormat format = params.getOutputFormat();
		DataSet data = pair.getFirstElement();
		int width = params.getWidth();
		double[] weights = pair.getSecondElement();
		
		boolean sortByWeights = true;
		double v = weights[0];
		for(int i=1;i<weights.length;i++){
			if(weights[i] != v){
				sortByWeights = false;
			}
		}
		
		int blockSpacer = params.getHeightOfSequenceLogo();
		int numBestForSorting = params.getNumberOfDependencies();
		if(numBestForSorting >=data.getElementLength()){
			numBestForSorting = data.getElementLength()-1;
		}
		
		
		
		prot.appendHeading( "Blocks of sequences" );
		
		int[] numPerChunk = params.getNumbersOfSequencesForBlocks();
		int[] blockHeight = params.getHeightsOfBlocks();
		int n = data.getNumberOfElements();
		for(int i=0;i<numPerChunk.length-1;i++){
			n -= numPerChunk[i];
		}

		if(n<=0){
			prot.appendWarning( "Block definition requires "+(-n)+" more sequences than present in the input data set. Terminating." );
			return null;
		}else{
			numPerChunk[numPerChunk.length-1] = n;

			for(int i=0;i<numPerChunk.length;i++){
				prot.append( "Block "+(i+1)+" contains "+numPerChunk[i]+" sequences and has height "+blockHeight[i]+".<br />" );
			}


			LinkedImageResult res = run(format, pdfPath, data, width, weights, numPerChunk, blockHeight, blockSpacer, numBestForSorting, sortByWeights );
			return res;
		}
	}
	
	
	private static LinkedImageResult run(OutputFormat format, String pdfPath, DataSet data, int width, double[] weights, int[] numPerChunk, int[] blockHeight, int blockSpacer, int numBestForSorting, boolean sortByWeights ) throws Exception{
		
		
		GraphicsAdaptor adaptor = GraphicsAdaptorFactory.getAdaptor( format );
	    
		plot( adaptor, data, width, weights.clone(), numPerChunk, blockHeight, blockSpacer, numBestForSorting, sortByWeights );
		
		//String pdfPath = ;
		
		pdfPath += "."+adaptor.getGraphicsExtension();//TODO
		
		adaptor.generateOutput( pdfPath );
		
		
		adaptor = GraphicsAdaptorFactory.getAdaptor( OutputFormat.PNG );
		
		blockSpacer = (int)Math.round(400.0/(double)width*blockSpacer);
		
		int[] blockHeightSmall = blockHeight.clone();
		for(int i=0;i<blockHeightSmall.length;i++){
			blockHeightSmall[i] = (int)Math.round(400.0/(double)width*blockHeightSmall[i]);
		}
		
		
		plot( adaptor, data, 600, weights.clone(), numPerChunk, blockHeightSmall, blockSpacer, numBestForSorting, sortByWeights );
		
		RasterizedAdaptor rast = (RasterizedAdaptor) adaptor;
		
		BufferedImage img = rast.getImage();
		
		LinkedImageResult res = new LinkedImageResult( "Dependency logo", "", img, new FileResult( format.name(), "Dependency logo for download", pdfPath ) );
		
		return res;
	}
	
	private static void plot(GraphicsAdaptor adaptor, DataSet data, int width, double[] weights, int[] numPerChunk, int[] blockHeight, int blockSpacer, int numBestForSorting, boolean sortByWeights) throws Exception {
		
		int height=SeqLogoPlotter.getHeightForDependencyLogo( data.getElementLength(), data.getNumberOfElements(), blockHeight, width, blockSpacer );

		Graphics2D g = adaptor.getGraphics( width, height );
		
		//Graphics2D g = (Graphics2D)img.getGraphics();
		g.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
		g.setRenderingHint(RenderingHints.KEY_RENDERING, RenderingHints.VALUE_RENDER_QUALITY);
		
		SeqLogoPlotter.plotDependencyLogo( data, null, 1, null, weights, g, width, 0, 0, numPerChunk, blockHeight, 0.03, blockSpacer, false, numBestForSorting, false, sortByWeights, true, 0.1);
	}
	
	
	

}
