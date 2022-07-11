package projects.talecorrect;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Date;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.TreeMap;
import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.DNADataSet;
import de.jstacs.data.alphabets.DiscreteAlphabet;
import de.jstacs.data.alphabets.GenericComplementableDiscreteAlphabet;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.data.sequences.annotation.SimpleSequenceAnnotationParser;
import de.jstacs.parameters.FileParameter;
import de.jstacs.parameters.Parameter;
import de.jstacs.parameters.FileParameter.FileRepresentation;
import de.jstacs.results.Result;
import de.jstacs.results.ResultSet;
import de.jstacs.results.TextResult;
import de.jstacs.tools.JstacsTool;
import de.jstacs.tools.ProgressUpdater;
import de.jstacs.tools.Protocol;
import de.jstacs.tools.ToolParameterSet;
import de.jstacs.tools.ToolResult;
import de.jstacs.tools.ui.cli.CLI;

public class CorrectTALESequences implements JstacsTool {
	public static void main( String[] args) throws Exception {
		CLI cli = new CLI(new CorrectTALESequences());
		cli.run(args);
	}
	
	@Override
	public ToolParameterSet getToolParameters() {
		
		LinkedList<Parameter> pars = new LinkedList<>();
		try {
			pars.add(new FileParameter("Sequences", "File with TALE-containing ONT assembly or extracted TALE sequences", "fasta,fa", true));
			pars.add(new FileParameter("N-terminus nHMMER File", "The output of nHMMER N-terminus", "txt", true));
			pars.add(new FileParameter("Repeats nHMMER File", "The output of nHMMER repeats", "txt", true));
			pars.add(new FileParameter("C-terminus nHMMER File", "The output of nHMMER C-terminus", "txt", true));
		}catch(Exception e){
			e.printStackTrace();
			return null;
		}
		
		return new ToolParameterSet(this.getShortName(), pars.toArray(new Parameter[0]));
	}
	@Override
	public ToolResult run(ToolParameterSet parameters, Protocol protocol, ProgressUpdater progress, int threads) throws Exception {		
		FileRepresentation nanoAssemblyWithTALEs = ((FileParameter)parameters.getParameterAt(0)).getFileContents();
		FileRepresentation nHmmerNFile = ((FileParameter)parameters.getParameterAt(1)).getFileContents();
		FileRepresentation nHmmerRFile = ((FileParameter)parameters.getParameterAt(2)).getFileContents();
		FileRepresentation nHmmerCFile = ((FileParameter)parameters.getParameterAt(3)).getFileContents();

		int countDifferneces_commonChar_vs_HomopolymerSubstitution=0;
		
		SimpleSequenceAnnotationParser parser = new SimpleSequenceAnnotationParser();
		System.out.println();
		DNADataSet ds = new DNADataSet( nanoAssemblyWithTALEs.getFilename(),'>', parser );
		
		String[] symbols_caseSensitive={"A", "C", "G", "T","a","c","g","t"};
		AlphabetContainer conCaseSensitive=new AlphabetContainer(new DiscreteAlphabet(false,symbols_caseSensitive));
		String[] symbols = {"A", "C", "G", "T", "-"};
		DiscreteAlphabet abc = new DiscreteAlphabet( true, symbols );

		int[] revComp = new int[symbols.length];
		revComp[0] = 3; //symbols[0]^rc = symbols[3]
		revComp[1] = 2; //symbols[1]^rc = symbols[2]
		revComp[2] = 1; //symbols[2]^rc = symbols[1]
		revComp[3] = 0; //symbols[3]^rc = symbols[0]
		revComp[4] = 4; //symbols[4]^rc = symbols[4]
		GenericComplementableDiscreteAlphabet abc2 = new GenericComplementableDiscreteAlphabet( true, symbols, revComp );
		AlphabetContainer con2=new AlphabetContainer( abc2 );
		HashMap<String, Sequence> seqTALEs=new HashMap<String, Sequence>();

		for(int i=0;i<ds.getNumberOfElements();i++){
			Sequence seq =ds.getElementAt(i);
			String header=seq.getAnnotation()[0].getResultAt(0).getValue().toString();
			String TALEname=header.split(" ")[0];
			seqTALEs.put(TALEname, seq);
		}
		ArrayList<nHMMERCorrection> correctionList=new ArrayList<nHMMERCorrection>();
		searchForCorrections(nHmmerNFile.getFilename(),0.75,200,seqTALEs,correctionList,3,con2,"N");
		searchForCorrections(nHmmerRFile.getFilename(),0.75,45,seqTALEs,correctionList,2,con2,"R");
		searchForCorrections(nHmmerCFile.getFilename(),0.75,200,seqTALEs,correctionList,3,con2,"C");
		
		System.out.println("candidate Corrections:");
		System.out.println(correctionList.toString());
		StringBuffer subList=new StringBuffer();
		subList.append("seqName"+"\t"+"position in uncorrected sequences"+"\t"+"type"+"\t"+"substitution"+"\n");
		for(String seqName: seqTALEs.keySet()){
			HashMap<Integer, Integer> pos_index=new HashMap<Integer,Integer>();
			int index=0;
			for(nHMMERCorrection nHcorr : correctionList){
				if(nHcorr.getSeqName().equals(seqName)){
					pos_index.put(nHcorr.getPosition(), index);
				}
				index++;
			}
			TreeMap<Integer, Integer> sortedPos_index = new TreeMap<>(pos_index);
			//System.out.println("Treemap: "+sortedPos_index);
			//The replacements/insertions are made from the end of the sequence to the beginning
			for(Integer pos:sortedPos_index.descendingKeySet()){
				System.out.println(correctionList.get(sortedPos_index.get(pos)).toString());
				int seqLength=seqTALEs.get(seqName).getLength();
				Sequence seqCorrected=null;
				System.out.println("pos: "+pos);
				if(correctionList.get(sortedPos_index.get(pos)).getType()=='i'){
					int correctionPosition=correctionList.get(sortedPos_index.get(pos)).getCorrectionPosition();
					
					String seqBefore=seqTALEs.get(seqName).getSubSequence(0,correctionPosition-1).toString();
					System.out.println("seq: "+seqTALEs.get(seqName).getSubSequence(correctionPosition-10,20).toString());
					//Homopolymer before pos
					String shortSeqBefore=seqTALEs.get(seqName).getSubSequence(correctionPosition-11,10).toString();
//					System.out.println("seqBefore: "+seqTALEs.get(seqName).getSubSequence(correctionPosition-11,10).toString());
					String homopolymerBefore=getHomopolymerBefore(shortSeqBefore.toUpperCase());
					//System.out.println(homopolymerBefore);
					String seqAfter= seqTALEs.get(seqName).getSubSequence(correctionPosition-1).toString();
					
					//homopolymer after pos
					String shortSeqAfter=seqTALEs.get(seqName).getSubSequence(correctionPosition-1,10).toString();
					String homopolymerAfter=getHomopolymerAfter(shortSeqAfter.toUpperCase());
//					System.out.println("homopolymerAfter: "+homopolymerAfter);
//					System.out.println("seqAfter"+": "+seqTALEs.get(seqName).getSubSequence(correctionPosition-1,10).toString());
					//char substitution=correctionList.get(sortedPos_index.get(pos)).getToChar();
					String longestHomopoylmer="";
					
					if(homopolymerAfter.charAt(0)==homopolymerBefore.charAt(0)){
						longestHomopoylmer=homopolymerAfter+homopolymerBefore;
					}
					else if(homopolymerAfter.length()>homopolymerBefore.length()){
						longestHomopoylmer=homopolymerAfter;
					}else if(homopolymerAfter.length()<homopolymerBefore.length()){
						longestHomopoylmer=homopolymerBefore;
					}else{
						longestHomopoylmer=homopolymerBefore;
						//System.err.println("Homopolymer before and after insertion are same length. Choose before.");
					}
					
					if(longestHomopoylmer.charAt(0)!=correctionList.get(sortedPos_index.get(pos)).getToChar()){
						System.err.println("Difference: HomopolymerChar ("+longestHomopoylmer.charAt(0)+") vs. CommonNucl ("+correctionList.get(sortedPos_index.get(pos)).getToChar()+")");
						
						//System.out.println("Difference: substitionHomopolymerChar vs substitionCommonNucl");
						countDifferneces_commonChar_vs_HomopolymerSubstitution++;
					}
					//subustitution for homoploymers with minimum length of 5 is homopolymer base, otherwise substition base is common base
					char substitution=correctionList.get(sortedPos_index.get(pos)).getToChar();
					if(longestHomopoylmer.length()>=5){
						substitution=longestHomopoylmer.charAt(0);
					}
					substitution=Character.toLowerCase(substitution);
					subList.append(seqName+"\t"+pos+"\t"+"insertion"+"\t"+correctionList.get(sortedPos_index.get(pos)).getFromChar()+" -> "+substitution+"\n");
					System.out.println("longestHomopoylmer: "+longestHomopoylmer);
					System.out.println("--seq Parts--");
					System.out.println("seqStart:");
					System.out.println(seqBefore.substring(seqBefore.length()-10));
					System.out.println("insertion");
					System.out.println(substitution);
					System.out.println("seqEnd");
					System.out.println(seqAfter.substring(0, 10));
					System.out.println();
					seqCorrected=Sequence.create(conCaseSensitive, seqBefore+substitution+seqAfter); 
				}else if(correctionList.get(sortedPos_index.get(pos)).getType()=='d'){
					subList.append(seqName+"\t"+pos+"\t"+"deletion"+"\t"+correctionList.get(sortedPos_index.get(pos)).getFromChar()+" -> "+"-"+"\n");
					String seqBefore=seqTALEs.get(seqName).getSubSequence(0,pos-1).toString();
					String seqAfter= seqTALEs.get(seqName).getSubSequence(pos).toString();
					seqCorrected=Sequence.create(conCaseSensitive, seqBefore+seqAfter); 
					System.out.println("deletion");
					System.out.println();
				}
				
				seqTALEs.replace(seqName, seqCorrected);
			}
		}
		
		
		StringBuffer correctedSequences =new StringBuffer();
		StringBuffer insertionList=new StringBuffer();
		StringBuffer igvToolsCountScript=new StringBuffer();
		igvToolsCountScript.append("#!/bin/bash\n");
		igvToolsCountScript.append("assembly=$1\n");
		igvToolsCountScript.append("bamFile=$2\n");
		igvToolsCountScript.append("pathFiles=$3\n");
		List<String> sortedKeys = new ArrayList<>(seqTALEs.keySet());
		Collections.sort(sortedKeys);
		for(String seqName : sortedKeys){
			//System.out.println(seqName);
			correctedSequences.append(">"+seqName+"\n");
			String seq=seqTALEs.get(seqName).toString();
			
			char[] acgt ={'a','c','g','t'};
			int pos=-1;
			int posBefore=-2;
			int countWigFiles=0;
			boolean isFirst=true;
			ArrayList<Integer> posList=new ArrayList<Integer>();
			for(char sub : acgt){
				for(int i=0;i<seq.length();i++){
					if((pos=seq.indexOf(sub, i)+1)>=1){
						//System.out.println(pos);
						insertionList.append(seqName+"\t"+pos+"\t"+sub+"\n");
						if((pos==posBefore+1)||(isFirst)){
							posList.add(pos);
							isFirst=false;
						}else{
							countWigFiles++;
							Collections.sort(posList);
							igvToolsCountScript.append("wigFile=$pathFiles"+"\"out.igvtools.count."+countWigFiles+".wig\"\n");
							//System.out.println(posList.toString()+"\t"+posList.size());
							//igvToolsCountScript.append("igvtools count -w 1 --bases --minMapQuality 1 --query=\""+seqName+":"+posList.get(0)+"-"+posList.get(posList.size()-1)+"\" $bamFile $wigFile $assembly &\n");
							igvToolsCountScript.append("igvtools count -w 1 --bases --minMapQuality 1 --query=\""+seqName+":"+posList.get(0)+"-"+posList.get(posList.size()-1)+"\" $bamFile $wigFile $assembly\n");
							posList.clear();
							posList.add(pos);
//							if(countWigFiles%5==0){
//								igvToolsCountScript.append("sleep(20)\n");
//							}
						}
						posBefore=pos;
						i=pos;
					}else{
						break;
					}
				}
			}
			
			
//			if(inputTyp.equals("TALEs")){
//				if(seq.subSequence(1, 4).equals("ATG")){
//					seq=seq.substring(1);
//				}
//				int startIndex=seq.indexOf("ATGGAT");
//				System.out.println("startIndex: "+startIndex);
//				if(startIndex==-1){
//					startIndex=0;
//					//warning
//					System.out.println("N-Terminus ATGGAT not found.. N-Terminus too short!!");
//				}
//				seq=seq.substring(startIndex);
//			}
			
			correctedSequences.append(seq+"\n");
		}
		System.out.println("countDifferneces commonChar vs HomopolymerSubstitution: "+countDifferneces_commonChar_vs_HomopolymerSubstitution);
	
		protocol.append( "\nWriting outputs.\n" );
		
		TextResult fr1 = new TextResult( "correctedTALEs", "Output with corrected sequences.", new FileRepresentation( "", correctedSequences.toString() ), "fa", getToolName(), "fasta/dna",true );
		TextResult fr2 = new TextResult( "insertionList_inCorrectedFile", "List of insertions with positions in corrected sequences.", new FileRepresentation( "", insertionList.toString() ), "tsv", getToolName(), "tsv",true );
		TextResult fr3 = new TextResult( "igvtools.count.Substitution.script", "Shell script to use for substition polishing.", new FileRepresentation( "", igvToolsCountScript.toString() ), "sh", getToolName(), "sh",true );
		TextResult fr4 = new TextResult( "substitionList", "List of substitutions with positions within uncorrected sequences.", new FileRepresentation( "", subList.toString() ), "tsv", getToolName(), "tsv",true );
		
		
		ResultSet set = new ResultSet( new Result[]{fr1,fr2,fr3,fr4} );
		return new ToolResult("Result of "+getToolName(), getToolName()+" on \""+nanoAssemblyWithTALEs.getFilename()+"\"", null, set, parameters, getToolName(), new Date(System.currentTimeMillis()) );
	}
	
	private static class nHMMERCorrection
	{
		private String seqName;
		private int position;
		private int correctionPosition;
		private char fromChar;
		private char toChar;
		private char type;
		
		public nHMMERCorrection( String seqName,  int position,  int correctionPosition,  char fromChar,  char toChar,  char type) {
			this(seqName, position, correctionPosition, fromChar, type);
			this.toChar = toChar;
		}
		
		public nHMMERCorrection( String seqName,  int position,  int correctionPosition,  char fromChar,  char type) {
			this.seqName = seqName;
			this.position = position;
			this.correctionPosition = correctionPosition;
			this.fromChar = fromChar;
			this.toChar = '-';
			this.type = type;
		}
		
		public String getSeqName() {
			return this.seqName;
		}
		
		public int getPosition() {
			return this.position;
		}
		
		public int getCorrectionPosition() {
			return this.correctionPosition;
		}
		
		public char getFromChar() {
			return this.fromChar;
		}
		
		public char getToChar() {
			return this.toChar;
		}
		
		@Override
		public String toString() {
			return this.seqName + " (" + this.correctionPosition + "): " + this.fromChar + " -> " + this.toChar;
		}
		
		public char getType() {
			return this.type;
		}
	}
	
	private static class nHmmerResult
	{
		private int startPosConsensus;
		private int startPos;
		private int endPos;
		private Sequence seq;
		private int strand;
		
		public nHmmerResult( Sequence seq,  int startPosConsensus,  int startPos,  int endPos,  int strand) {
			this.strand = 0;
			this.seq = seq;
			this.startPosConsensus = startPosConsensus;
			this.startPos = startPos;
			this.endPos = endPos;
			if (endPos - startPos < 0) {
				this.strand = -1;
			}
			else {
				this.strand = 1;
			}
		}
		
		public int getStartPosConsensus() {
			return this.startPosConsensus;
		}
		
		public int getStartPos() {
			return this.startPos;
		}
		
		public int getEndPos() {
			return this.endPos;
		}
		
		public Sequence getSequence() {
			return this.seq;
		}
		
		public int getStrand() {
			return this.strand;
		}
		
		@Override
		public String toString() {
			 String result = this.startPos + " " + this.seq + " " + this.endPos;
			return result;
		}
	}
	
	private static String skipLinesBufferedReader(BufferedReader BR, int numLines) throws IOException{
		String line="";
		for(int s=0;s<numLines;s++){
			line = BR.readLine();
		}
		return line;
	}
	
	private static String getHomopolymerBefore (String SeqBeforeInsertion){
		char charAtEnd= SeqBeforeInsertion.charAt(SeqBeforeInsertion.length()-1);
		String homopolymer=SeqBeforeInsertion.substring(SeqBeforeInsertion.length()-1, SeqBeforeInsertion.length());
		int i=SeqBeforeInsertion.length()-1;
		while(SeqBeforeInsertion.charAt(i)==charAtEnd){
			//System.out.println(i+": "+SeqBeforeInsertion.charAt(i));
			homopolymer=SeqBeforeInsertion.substring(i,SeqBeforeInsertion.length());
			i--;
		}
			return homopolymer;
	}
	
	private static String getHomopolymerAfter (String SeqAfterInsertion){
		char charAtStart= SeqAfterInsertion.charAt(0);
		String homopolymer=SeqAfterInsertion.substring(0, 1);
		int i=0;
		while((SeqAfterInsertion.charAt(i)==charAtStart)&&(i<SeqAfterInsertion.length())){
			//System.out.println(i+": "+SeqAfterInsertion.charAt(i));
			homopolymer=SeqAfterInsertion.substring(0,i+1);
			i++;
		}
			return homopolymer;
	}

	private static void searchForCorrections(String nHmmerPath,double minAcc,int minScore,HashMap<String, Sequence> seqTALEs,ArrayList<nHMMERCorrection>correctionList,int maxConitgousCorrection, AlphabetContainer con2, String type) throws Exception{
		BufferedReader BR = new BufferedReader(new FileReader(nHmmerPath));
		
		String line=null;
		boolean readHeader=true;
		boolean parseResult=false;
		String aktResultConsensusSeq="";
		int aktResultStartPosConsensusSeq=-1;
		String aktResultPosteriorProbability="";
		String nHmmerOutSeqname="";

		while ((line = BR.readLine()) != null){
			
			if(readHeader){
				if(line.startsWith("Annotation for each hit  (and alignments):")){
					readHeader=false;
				}
			}else{
				
				if(line.equals("")){
					line=skipLinesBufferedReader(BR,1);
				}
				int hmmTo=-1;
				if(line.startsWith(">>")){
					parseResult=false;
					nHmmerOutSeqname=line.split(" ")[1];
					//System.out.println(nHmmerOutSeqname);
					for(int s=0;s<3;s++){
						line = BR.readLine();
					}
					line=line.trim().replaceAll(" +"," ");
					String[] lineSplit=line.split(" ");
					hmmTo=Integer.parseInt(lineSplit[5]);

					double acc=Double.parseDouble(lineSplit[14]);
					double score=Double.parseDouble(lineSplit[1]);
					if((acc>=minAcc)&&(score>=minScore)){ 
						parseResult=true;
						line=skipLinesBufferedReader(BR,3);
					}
				}else{
					if(parseResult){
						while(line.equals("")){
							line=skipLinesBufferedReader(BR,1);
						}
						if(line.startsWith("Internal pipeline statistics summary")){
							break;
						}
						
						line=line.trim();
						line=line.replaceAll(" +"," ");
							
						aktResultConsensusSeq=line.split(" ")[2];
		
						aktResultStartPosConsensusSeq=Integer.parseInt(line.split(" ")[1]);
						//System.out.println(aktResultConsensusSeq);
						line=skipLinesBufferedReader(BR,2);
						line=line.trim().replaceAll(" +"," ");
						
						//System.out.println("line: "+line);
							if(!line.split(" ")[1].equals("-")){
								int startPosTALE=Integer.parseInt(line.split(" ")[1]);
								int endPosTALE=Integer.parseInt(line.split(" ")[3]);
								int strand=0;
								if(endPosTALE-startPosTALE<0){
									strand=-1; //match on "- strand"
								}else{
									strand=1; //match on "+ strand"
								}
								
								Sequence subSeqTALEnHmmer=Sequence.create(con2, line.split(" ")[2]);
								
								nHmmerResult aktResult=new nHmmerResult(subSeqTALEnHmmer,aktResultStartPosConsensusSeq, startPosTALE, endPosTALE,strand);
								
								//read line with Posterior Probability
								line = BR.readLine();
								line=line.trim().replaceAll(" +"," ");
								
								aktResultPosteriorProbability=line;
								
								if(aktResultPosteriorProbability.contains(".")){
//									System.out.println(aktResultConsensusSeq);
//									System.out.println(aktResult.getSequence());
//									System.out.println(aktResultPosteriorProbability);
									ArrayList<Integer> indexDifferences= new ArrayList<Integer>();
									
									for(int p=aktResultPosteriorProbability.indexOf(".");p<aktResultPosteriorProbability.length();p++){
									
										if(((!((aktResult.getStartPosConsensus()+p)>(864)))&&(type=="N"))||(type!="N")){ //consensus seq N-terminus not correct the last 10 bases, as these belong to 1. repeat and are corrected separately
											if((aktResultPosteriorProbability.charAt(p)=='.')&&(aktResult.getSequence().toString().charAt(p)=='-')){
												indexDifferences.add(p);
											}
											if((aktResultPosteriorProbability.charAt(p)=='.')&&(aktResult.getSequence().toString().charAt(p)!='-')){
												//throw new Exception("Correction type is not an insertion!");
											}
										}
										
									}
									//System.out.println("PotentialToCorrectPositions: "+indexDifferences.toString());
									
									int countContiguousPositions=1;
									for(int p=0;p<indexDifferences.size()-1;p++){
										int indexP=indexDifferences.get(p);
										// check for long contiguous corrections
										if((indexP+1)==indexDifferences.get(p+1)){
											
											countContiguousPositions++;
											
										}
										if(((indexP+1)!=indexDifferences.get(p+1))||(p==indexDifferences.size()-2)){
											//delete too long contiguous streches!
											if(countContiguousPositions>maxConitgousCorrection){
												for(int d=0;d<countContiguousPositions;d++){
													indexDifferences.remove(p+1-d);
												}
												
											}
											countContiguousPositions=1;
										}
									}
									//System.out.println("PotentialToCorrectPositions: "+indexDifferences.toString());
									
									//System.out.println("Hier:"+nHmmerOutTALEname);
									//correction in fasta File!
									//System.out.println(aktResult.getStartPos());
									
									for(int p=0;p<indexDifferences.size();p++){
										int posInNHMMER=indexDifferences.get(p);
//										System.out.println(posInNHMMER);//counted from 0
//										System.out.println(aktResult.getSequence().toString().substring(posInNHMMER, posInNHMMER+1));
										int posInInputSeq;
										int posInInputSeqWithP;
										if(aktResult.getStrand()==-1){
											posInInputSeq=aktResult.getStartPos()+((posInNHMMER-p)*aktResult.getStrand()+1);
											posInInputSeqWithP=aktResult.getStartPos()+((posInNHMMER)*aktResult.getStrand()+1);
										}else{
											posInInputSeq=aktResult.getStartPos()+((posInNHMMER-p)*aktResult.getStrand());
											posInInputSeqWithP=aktResult.getStartPos()+((posInNHMMER)*aktResult.getStrand());
										}
										
//										System.out.println(posInInputSeq);
//										System.out.println(seqTALEs.get(nHmmerOutSeqname).getSubSequence(posInInputSeq, 1));
										//int seqLength=seqTALEs.get(nHmmerOutTALEname).getLength();
										char fromChar='N';
										String substitution=aktResultConsensusSeq.substring(posInNHMMER, posInNHMMER+1).toUpperCase();
										Sequence substition_seq=null;
										if((aktResult.getStrand())==-1){
											fromChar=aktResult.getSequence().complement().toString().substring(posInNHMMER, posInNHMMER+1).charAt(0);
											substition_seq=Sequence.create(con2, substitution).complement();
										}else{
											fromChar=aktResult.getSequence().toString().substring(posInNHMMER, posInNHMMER+1).charAt(0);
											substition_seq=Sequence.create(con2, substitution);
										}
										 
										
//										System.out.println(nHmmerOutSeqname);
//										System.out.println("FromChar:"+fromChar);
//										System.out.println("ToChar:"+substitution.charAt(0));
										correctionList.add(new nHMMERCorrection(nHmmerOutSeqname,posInInputSeqWithP,posInInputSeq,fromChar,substition_seq.toString().charAt(0),'i'));			
								
									}
									
								}
								//Deletions
								if(aktResultConsensusSeq.matches("[acgt]*[.][acgt]*")||aktResultConsensusSeq.matches("[acgt]*[..][acgt]*")){
//									System.out.println(aktResultConsensusSeq);
//									System.out.println(aktResult.getSequence());
//									System.out.println(aktResultPosteriorProbability);
									ArrayList<Integer> indexDifferences= new ArrayList<Integer>();
									
									
									for(int p=aktResultConsensusSeq.indexOf(".");p<aktResultConsensusSeq.length();p++){
										if(aktResultConsensusSeq.charAt(p)=='.'){
											indexDifferences.add(p);
										}
									}
									for(int p=0;p<indexDifferences.size();p++){
										int posInNHMMER=indexDifferences.get(p);
//										System.out.println(posInNHMMER);//counted from 0
//										System.out.println(aktResult.getSequence().toString().substring(posInNHMMER, posInNHMMER+1));
										int posInInputSeq=aktResult.getStartPos()+(posInNHMMER*aktResult.getStrand());

//										System.out.println(posInInputSeq);
//										System.out.println(seqTALEs.get(nHmmerOutSeqname).getSubSequence(posInInputSeq, 1));
										//int seqLength=seqTALEs.get(nHmmerOutTALEname).getLength();
									
										char fromChar='N';
										if((aktResult.getStrand())==-1){
											fromChar=aktResult.getSequence().complement().toString().substring(posInNHMMER, posInNHMMER+1).charAt(0);
										}else{
											fromChar=aktResult.getSequence().toString().substring(posInNHMMER, posInNHMMER+1).charAt(0);
										}
//										System.out.println(nHmmerOutSeqname);
//										System.out.println("FromChar:"+fromChar);

										correctionList.add(new nHMMERCorrection(nHmmerOutSeqname,posInInputSeq,posInInputSeq,fromChar,'d'));			
								
									}
									
									//System.err.println("WARNING: Correction type is a deletion! Please check: "+nHmmerOutSeqname+": "+aktResult.toString());
								}
							}else{
								line=skipLinesBufferedReader(BR,1);
							}
							
						
						
						//parseResult=false;
						
					}
				}
			}
			
		}
		BR.close();
	}
	
	@Override
	public String getToolName() {
		return "CorrectTALESequences";
	}
	@Override
	public String getToolVersion() {
		return "0.1";
	}
	@Override
	public String getShortName() {
		return "correct";
	}
	@Override
	public String getDescription() {
		return "Corrects TALE sequences.";
	}
	@Override
	public String getHelpText() {
		// TODO Auto-generated method stub
		return null;
	}
	@Override
	public ResultEntry[] getDefaultResultInfos() {
		// TODO Auto-generated method stub
		return null;
	}
	@Override
	public ToolResult[] getTestCases(String path) {
		// TODO Auto-generated method stub
		return null;
	}
	@Override
	public void clear() {
		// TODO Auto-generated method stub
		
	}
	@Override
	public String[] getReferences() {
		// TODO Auto-generated method stub
		return null;
	}

}




