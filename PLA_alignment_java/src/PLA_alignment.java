/**
 * Java program for PLA alignment from raw reads
 * Command line argument example: input="abc\def.txt" (no whitespaces are allowed around the equal signs)
 */


// Import packages
import java.io.*;
import java.util.Arrays;
import java.util.ArrayList;
import java.util.List;
import java.util.Set;
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.TreeSet;
import java.util.Objects;
import java.util.zip.*;

//import java.lang.Math;
import java.time.format.DateTimeFormatter;
import java.time.LocalDateTime;

import org.apache.commons.lang3.StringUtils;
import org.apache.commons.text.similarity.LevenshteinDistance;
import com.google.common.collect.*;


public class PLA_alignment
{
	
	// Calculate Hamming distance
	// Can allow N, which matches to all bases
	public static int HammingDistanceCalculator(String str1, String str2, boolean allowN)
	{
		if (str1.length() != str2.length()) {throw new IllegalArgumentException("Input strings must have equal lengths!");}
		else
		{
			int dist = 0;
			
			for (int i=0; i<str1.length(); i++)
			{
				
				if (str1.charAt(i) != str2.charAt(i) )
				{
					if (allowN && ((str1.charAt(i)=='N') || (str2.charAt(i)=='N'))) {continue;} else {dist++;}
				}
			}
			return dist;
		}
	}
	public static int HammingDistanceCalculator(String str1, String str2) {return HammingDistanceCalculator(str1, str2, false);} // overloading: default allowN is false
	
	public static void main(String[] args)
	{
		// Time format
		DateTimeFormatter time_formatter = DateTimeFormatter.ofPattern("yyyy-MM-dd HH:mm:ss");
		
		switch (args[0])
		{
		
		/**
		 * Trim the first N bases from a fastq.gz file
		 * I=... O=... N=...
		 * This is used to remove the mosaic sequence from the RNA library if RNA and PLA are sequenced together
		 * A read will only be trimmed if it is longer than N bases
		 * 
		 * Input arguments:
		 * 		I: path to input file (fastq.gz format)
		 * 		O: path to store output file (fastq.gz format)
		 * 		N: number of bases to trim
		 */
		case "ReadTrimming":
		{
			// Parse the arguments
			String I = "", O = "";
			int N = 0;
			for (int i=1; i<args.length; i++)
			{
				String[] j = args[i].split("=");
				
				if (j.length < 2) {throw new IllegalArgumentException("ReadTrimming: Whitespaces are not allowed between argument key specifier and the argument!");}
				
				switch (j[0])
				{
				case "I": I = j[1]; break;
				case "O": O = j[1]; break;
				case "N": N = Integer.parseInt(j[1]); break;
				default: throw new IllegalArgumentException("ReadTrimming: Invalid argument key specifier!");
				}
			}
			
			try (
					BufferedReader brI = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(I)))); // input file
					BufferedWriter bwout = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(O)))) // output file
				 )
			{
				// Write out the command line arguments
				System.out.println();
				for (String j : args)
				{
					System.out.println(j);
				}
				System.out.println();
				
				// Initialize variables
				String line;
				int counter = 0; // line number counter
				int short_counter = 0; // counter for number of reads that are shorter than N
				long my_timer = System.currentTimeMillis(); // timer
				System.out.printf("%s   ReadTrimming   Start%n", LocalDateTime.now().format(time_formatter));
				while ((line=brI.readLine()) != null) 
				{
					if (((counter % 4) == 1) || ((counter % 4) == 3)) // trim sequence and quality score lines
					{
						if (line.length() > N)
						{
							bwout.write(line.substring(N));
							bwout.newLine();
						}
						else
						{
							bwout.write(line);
							bwout.newLine();
							short_counter++;
						}
						
						if (((counter % 4) == 1) && ((((counter-1)/4+1) % 1000000) == 0))
						{
							System.out.printf("%s   ReadTrimming   Processed %,15d reads   Elapsed time for last 1,000,000 reads: %ds%n",
									LocalDateTime.now().format(time_formatter), (counter-1)/4+1, (System.currentTimeMillis()-my_timer)/1000);
							my_timer = System.currentTimeMillis();
						}
					}
					else
					{
						bwout.write(line);
						bwout.newLine();
					}
					
					counter++;
				}
				System.out.printf("%s   ReadTrimming   Done: processed %,d reads%n", LocalDateTime.now().format(time_formatter), (counter-1)/4+1);
				System.out.printf("\tThere are %,d reads shorter than %d%n", short_counter, N);
				
			} catch (IOException e) { throw new IllegalArgumentException("Invalid file paths!");}
			break;
			
		}
		
		/**
		 * Alignment from raw reads for drop-seq runs
		 * R1=... R2=... O=... ABfile=... SUMMARY=... HEADER=... (for doing alignment of paired reads 1 by 1)
		 * 
		 * Input arguments:
		 * 		R1: path to read 1 file (fastq.gz format)
		 * 		R2: path to read 2 file (fastq.gz format)
		 * 		ABfile: path to PLA target-DNA barcode lookup table (csv format)
		 * 		O: path to store output file (txt.gz format)
		 * 		SUMMARY: directory to store summary files (txt format) (default is current working directory)
		 * 		HEADER: whether the ABfile has header to be skipped (default is false)
		 * 
		 * Output format: cell barcodes , UMI , AB1 ID , AB2 ID
		 */
		case "ReadAlignmentDropSeq":
		{
			
			// Parse the arguments
			String R1 = "", R2 = "", ABfile = "", O = "";
			String SUMMARY = System.getProperty("user.dir") + File.separator + "ReadAlignmentDropSeq_summary.txt"; // default summary file directory
			boolean skip_header=false;
			for (int i=1; i<args.length; i++)
			{
				String[] j = args[i].split("=");
				
				if (j.length < 2) {throw new IllegalArgumentException("Whitespaces are not allowed between argument key specifier and the argument!");}
				
				switch (j[0])
				{
				case "R1": R1 = j[1]; break;
				case "R2": R2 = j[1]; break;
				case "ABfile": ABfile = j[1]; break;
				case "O": O = j[1]; break;
				case "SUMMARY": SUMMARY = j[1]; break;
				case "HEADER": skip_header = "true".equalsIgnoreCase(j[1]); break;
				default: throw new IllegalArgumentException("ReadAlignmentDropSeq: Invalid argument key specifier!");
				}
			}
		
			try (
					BufferedReader br1 = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(R1)))); // read 1
					BufferedReader br2 = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(R2)))); // read 2
					BufferedReader brAB = new BufferedReader(new FileReader(ABfile)); // AB-DNA barcode look up table
					BufferedWriter bwout = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(O)))); // output file
					BufferedWriter bwsum = new BufferedWriter(new FileWriter(SUMMARY)) // summary file
				 )
			{
				
				// Write out the command line arguments
				System.out.println();
				for (String j : args)
				{
					bwsum.write(j); bwsum.newLine();
					System.out.println(j);
				}
				bwsum.newLine();
				System.out.println();
				
				// Read the AB look up table into an array
				List<List<String>> ABarray = new ArrayList<List<String>>();
				String ABline;
				if (skip_header) {brAB.readLine();} // skip first line
				while ((ABline=brAB.readLine()) != null)
				{
					String[] values = ABline.split(",");
					ABarray.add(Arrays.asList(values));
				}
	
				
				/**
				 * Process read 1 and 2
				 * Read 1 contains cell barcode and UMI
				 * Read 2 contains PLA products
				 * Output format: <cell barcode> , <UMI> , <AB1 ID> , <AB2 ID>
				 */
				
				// Set up the counters for the summary files
				int short_read_counter = 0; // counter for the number of reads with fewer than 75 bases
				int badUMI_counter = 0; // counter for the number of reads with invalid UMIs due to repeated G
				int excessiveG_counter = 0; // counter for the number of reads with too many Ns
				int excessiveN_counter = 0; // counter for the number of reads with too many Ns
				int bad_connector_counter = 0; // counter for the number of reads with non-matching connector sequence
				int non_matching_AB_counter = 0; // counter of the number of reads with non-matching AB barcode
				
				// Set up for alignment
				String line1, line2;
				int counter = 0; // line number counter
				int PLA_counter = 0; // PLA product counter
				String connector = "TCGTGTCGTGTCGTGTCTAAAG"; // connector sequence
				
				// Start reading
				long my_timer = System.currentTimeMillis();
				System.out.printf("%s   ReadAlignmentDropSeq   Start alignment%n", LocalDateTime.now().format(time_formatter));
				bwsum.write("Start alignment at " + LocalDateTime.now()); bwsum.newLine();
				
				while (((line1=br1.readLine()) != null) && ((line2=br2.readLine()) != null))
				{
					if ((counter % 4) == 0)
					{
						// Make sure that the identifiers of read 1 and 2 match up
						if (!Objects.equals(line1.split(" ")[0], line2.split(" ")[0]))
						{
							throw new IOException("Read 1 (" + String.format("%s",line1.split(" ")[0]) + ") and read 2 (" + String.format("%s",line2.split(" ")[0]) + ") do not match.");
						}
					}
							
					else if ((counter % 4) == 1)
					{
						
						
						
						// Expected location of connector starts at index 39 (base 40th)
						int connector_start = 39;					
						int connector_start_temp = 39; // temporary starting location of connector, for used in for loop
						
						
						// Skip reads with at least 7 occurrences of G in the UMI region
						if (line2.length() < 75)
						{
							short_read_counter++;
						}
						
						else if (StringUtils.countMatches(line1.substring(12), "G") >= 7)
						{
							badUMI_counter++;
						}
						
						// Skip reads with excessive Gs (reads that contain only G in the connector region)
						else if ( ( StringUtils.countMatches(line2.substring(connector_start,connector_start+connector.length()), "G") ) == connector.length() )
						{
							excessiveG_counter++;
						}
						
						// Check if read 2 has more than 10 N's
						else if (StringUtils.countMatches(line2, "N") > 10)
						{
							excessiveN_counter++;
						}
						
						else
						{							
							// Locate the connector using Levenshtein distance, max allowed distance is 3
							// Does not need to take into account N, since the connector region doesn't usually contain N
							int[] connector_shift = new int[] {0, -1, 1}; // allows 1-base shift
							boolean match_connector = false; // true if matching connector is found
							
							int connector_distance = 2; // lowest Levenshtein distance found between the true connector and the test connector sequence
							for (int shift_i : connector_shift)
							{
								int temp_distance = LevenshteinDistance.getDefaultInstance().apply(line2.substring(connector_start+shift_i, connector_start+shift_i+connector.length()), connector);
								if ((temp_distance <= 2) && (temp_distance < connector_distance))
								{
									connector_distance = temp_distance;
									connector_start_temp = connector_start + shift_i;
									match_connector = true;
									
									if (temp_distance == 0)
									{ break; }
								}
							}
							connector_start = connector_start_temp;
							
							if (match_connector)
							{
								
								// Initialize the AB ID
								String AB1_ID = "Unknown";
								String AB2_ID = "Unknown";
							
								// Initialize variables to store hamming distance for each AB1 and 2 barcode
								int[] temp1 = new int[ABarray.size()];
								int[] temp2 = new int[ABarray.size()];
								
								// Check if there is a frameshift, in order to locate AB2 correctly
								int shift_j = line2.substring(connector_start, connector_start+connector.length()+1).indexOf("TAAAG"); // location of AAAG in the found connector region
								if (shift_j == -1)
								{
									if ((((counter-1)/4+1) % 1000000) == 0)
									{
										System.out.printf("%s   ReadAlignmentDropSeq   Processed %,15d records   Elapsed time for last 1,000,000 reads: %ds%n",
												LocalDateTime.now().format(time_formatter), (counter-1)/4+1, (System.currentTimeMillis()-my_timer)/1000);
										my_timer = System.currentTimeMillis();
									}
									
									counter++;
									bad_connector_counter++;
									continue;
								}
								shift_j = 17 - shift_j; // number of bases to shift to the left
								
								// Found AB barcodes
								String AB1_found = line2.substring(connector_start-18, connector_start-18+8);
								String AB2_found = line2.substring(connector_start+25-shift_j, connector_start+25-shift_j+8);
								
								// Calculate Hamming distance only if there is at most 1 N in the found barcode
								if ((StringUtils.countMatches(AB1_found, "N") <= 1) && (StringUtils.countMatches(AB2_found, "N") <= 1))
								{
									// Counter for # of matches with 1 hamming distance
									int match_counter1 = 0, match_counter2 = 0;
									// Index for matches with 1 hamming distance
									int match_index1 = -1, match_index2 = -1;
									
									for (int i=0; (i<ABarray.size()) && ((AB1_ID=="Unknown") || (AB2_ID=="Unknown")) && ((match_counter1<2) && (match_counter2<2)); i++)
									{
										// Calculate Hamming distance
										temp1[i] = HammingDistanceCalculator(AB1_found, ABarray.get(i).get(1), true);
										temp2[i] = HammingDistanceCalculator(AB2_found, ABarray.get(i).get(1), true);
										
										// Allow early termination of for loop if found an exact match
										if (temp1[i] == 0)
										{
											AB1_ID = ABarray.get(i).get(0);
										}
										else if (temp1[i] == 1) {match_counter1++; match_index1 = i;}
										
										if (temp2[i] == 0)
										{
											AB2_ID = ABarray.get(i).get(0);
										}
										else if (temp2[i] == 1) {match_counter2++; match_index2 = i;}
										
									}
									
									// Find unambiguous match with 1 Hamming distance (ie, discard reads that have more than 1 matches with 1 hamming distance
									if (match_counter1 == 1) {AB1_ID = ABarray.get(match_index1).get(0);}
									if (match_counter2 == 1) {AB2_ID = ABarray.get(match_index2).get(0);}
									
									if (!Objects.equals(AB1_ID,"Unknown") && !Objects.equals(AB2_ID,"Unknown"))
									{
										bwout.write(line1.substring(0, 12)+","+line1.substring(12)+","+AB1_ID+","+AB2_ID);
										bwout.newLine();
										PLA_counter++;
									}
									else
									{
										non_matching_AB_counter++;
									}
								}
							}
							
							else
							{
								bad_connector_counter++;
							}
							
							
							
						}
						
						
						
						
						if ((((counter-1)/4+1) % 1000000) == 0)
						{
							System.out.printf("%s   ReadAlignmentDropSeq   Processed %,15d records   Elapsed time for last 1,000,000 reads: %ds%n",
									LocalDateTime.now().format(time_formatter), (counter-1)/4+1, (System.currentTimeMillis()-my_timer)/1000);
							my_timer = System.currentTimeMillis();
						}
					}
					
					
//					if ((((counter-1)/4)+1)>2000000) {System.out.printf("%s   %s   Processed %,d lines%n", LocalDateTime.now().format(time_formatter), counter); break;} // for testing purposes
				
					counter++;
				}
				
				
				System.out.printf("%s   ReadAlignmentDropSeq   Done: processed %,d reads%n", LocalDateTime.now().format(time_formatter), (counter-1)/4+1);
				System.out.printf("\tNumber of valid PLA products: %,15d%n", PLA_counter);
				System.out.println();
				
				// Write to summary file
				bwsum.write("ReadAlignmentDropseq: Finished at " + LocalDateTime.now().withNano(0) + ", processed " + String.format("%,d",(counter-1)/4+1) + " reads"); bwsum.newLine();
				bwsum.write("Number of valid PLA products: " + String.format("%,d", PLA_counter)); bwsum.newLine();
				bwsum.write("Number of records discarded because of read 2 being too short: " + String.format("%,d",short_read_counter)); bwsum.newLine();
				bwsum.write("Number of records discarded because of excessive G in UMIs: " + String.format("%,d",badUMI_counter)); bwsum.newLine();
				bwsum.write("Number of records discarded because of excessive number of Ns: " + String.format("%,d",excessiveN_counter)); bwsum.newLine();
				bwsum.write("Number of records discarded because of excessive G in read 2: " + String.format("%,d",excessiveG_counter)); bwsum.newLine();
				bwsum.write("Number of records discarded because of non-matching connector sequence: " + String.format("%,d",bad_connector_counter)); bwsum.newLine();
				bwsum.write("Number of records discarded because of non-matching antibody barcode: " + String.format("%,d",non_matching_AB_counter)); bwsum.newLine();

			} catch (IOException e) { throw new IllegalArgumentException("Invalid file paths!");}
			break;
			
		}
			
		/**
		 * Alignment from raw reads for smart-seq runs
		 * R1List=... O=... ABfile=... SUMMARY=... HEADER=...
		 * 
		 * Input arguments:
		 * 		R1List: path to a csv file containing the read 1 files (fastq.gz format) and the corresponding cell barcode
		 * 		^^^^^^ path/to/R1.fastq.gz,cell-barcode
		 * 		O: path to store output file (txt.gz format)
		 * 		ABfile: path to PLA target-DNA barcode lookup table (csv format)
		 * 		SUMMARY: directory to store summary files (txt format) (default is current working directory)
		 * 		HEADER: whether the ABfile has header to be skipped (default is false)
		 * 
		 * Output format: cell barcodes,UMI,AB1 ID,AB2 ID
		 * 
		 * The summary file also contains the found AB barcode and their read counts (only keep those with read counts >= 100)
		 */
		case "ReadAlignmentSmartSeq":
		{
			
			// Parse the arguments
			String R1List = "", ABfile = "", O = "";
			String SUMMARY = System.getProperty("user.dir") + File.separator + "ReadAlignmentSmartSeq_summary.txt"; // default summary file directory
			boolean skip_header = false;
			for (int i=1; i<args.length; i++)
			{
				String[] j = args[i].split("=");
				
				if (j.length < 2) {throw new java.lang.IllegalArgumentException("Whitespaces are not allowed between argument key specifier and the argument!");}
				
				switch (j[0])
				{
				case "R1List": R1List = j[1]; break;
				case "O": O = j[1]; break;
				case "ABfile": ABfile = j[1]; break;
				case "SUMMARY": SUMMARY = j[1]; break;
				case "HEADER": skip_header = "true".equalsIgnoreCase(j[1]); break;
				default: throw new java.lang.IllegalArgumentException("ReadAlignmentSmartSeq: Invalid argument key specifier!");
				}
			}
			
//			System.out.println(R1List);
		
			try (
					BufferedReader br1List = new BufferedReader(new FileReader(R1List)); // List of input read 1 files and the cell barcodes 
					BufferedReader brAB = new BufferedReader(new FileReader(ABfile)); // AB-DNA barcode look up table
					BufferedWriter bwout = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(O)))); // output file
					BufferedWriter bwsum = new BufferedWriter(new FileWriter(SUMMARY)) // summary file
				 )
			{
				
				// Write out the command line arguments
				System.out.println();
				for (String j : args)
				{
					bwsum.write(j); bwsum.newLine();
					System.out.println(j);
				}
				bwsum.newLine();
				System.out.println();
			
				
				// Read the AB look up table into an array
				List<List<String>> ABarray = new ArrayList<List<String>>();
				String ABline;
				if (skip_header) {brAB.readLine();} // skip first line
				while ((ABline=brAB.readLine()) != null)
				{
					String[] values = ABline.split(",");
					ABarray.add(Arrays.asList(values));
				}
	
				// Hash Multiset to store unique UMIs and their number of occurrences
				Multiset<String> ABcounts = HashMultiset.create();
				
				/**
				 * Process read 1
				 * Read 1 contains PLA products
				 * UMI region: 2nd base to 17th base (16-base long)
				 * Output format: <cell barcode> , <UMI> , <AB1 ID> , <AB2 ID>
				 * Cell barcode is equal to sample barcode
				 */

				// Set up the (total) counters for all read 1 files for the summary file
				
				int read_counter = 0; // read number counter for all fastq files
				int PLA_counter = 0; // PLA product counter
				int bad_connector_counter = 0; // counter for reads with non-matching connector
				int bad_UMI_counter = 0; // counter of reads with bad UMI region
				int non_matching_AB_counter = 0; // counter of the number of reads with non-matching AB barcode
				
				// Set up for alignment
				String R1List_temp; // current read 1 file
				String connector = "TCGTGTCGTGTCGTGTCTAAAG"; // connector sequence
				String[] R1List_temp_array; // array to store each line in R1List
				
				// Start reading
				long my_timer = System.currentTimeMillis();
				System.out.printf("%s   Start alignment%n", LocalDateTime.now().format(time_formatter));
				bwsum.write("Start alignment at " + LocalDateTime.now()); bwsum.newLine(); bwsum.newLine();
				
				while ((R1List_temp=br1List.readLine()) != null)
				{
					R1List_temp_array = R1List_temp.split(",");
					
					
					try (
							BufferedReader br1 = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(R1List_temp_array[0])))); // read1
						)
					{
						// Set up the counters for each read 1 file for the summary file
						int counter = 0; // line number counter in fastq file
//						int PLA_counter_temp = 0; // PLA product counter
//						int bad_connector_counter_temp = 0; // counter for reads with non-matching connector
//						int bad_UMI_counter_temp = 0; // counter of reads with bad UMI region
//						int non_matching_AB_counter_temp = 0; // counter of the number of reads with non-matching AB barcode
						
						String line1; // current line in current read 1 file
						while ((line1=br1.readLine()) != null)
						{
							if ((counter % 4) == 1)
							{			
								
								read_counter++;
								
								// Check if the UMI region has an excessive amount of A
								if (StringUtils.countMatches(line1.substring(1, 17), "A") >= 10)
								{
									if ((read_counter % 1_000_000) == 0)
									{
										System.out.printf("%s   ReadAlignmentSmartSeq   Processed %,15d reads   Elapsed time for last 1,000,000 reads: %ds%n",
												LocalDateTime.now().format(time_formatter), read_counter, (System.currentTimeMillis()-my_timer)/1000);
										my_timer = System.currentTimeMillis();
									}
									
									counter++;
									bad_UMI_counter++;
//									bad_UMI_counter_temp++;
									continue;
								}
								
								int connector_start = 39; // expected starting location of connector is at index 39 (base 40th)
								int connector_start_temp = 39; // temporary starting location of connector, for used in for loop
								
								// Locate the connector using Levenshtein distance, max allowed distance is 2
								// Does not need to take into account N, since the connector region doesn't usually contain N
								int[] connector_shift = new int[] {0, -1, 1}; // only allow up to 1-base shift of the expected starting location
								boolean match_connector = false; // true if matching connector is found
								int connector_distance = 3; // lowest Levenshtein distance found between the true connector and the test connector sequence, will be updated during the for loop
								for (int shift_i : connector_shift)
								{
									int temp_distance = LevenshteinDistance.getDefaultInstance().apply(line1.substring(connector_start+shift_i, connector_start+shift_i+connector.length()), connector);
									if ((temp_distance <= 2) && (temp_distance < connector_distance))
										// match within 2 Levenshtein distance, and only if the temp_distance is lower than previous iterations
									{
										connector_distance = temp_distance;
										connector_start_temp = connector_start + shift_i;
										match_connector = true;
										
										if (temp_distance == 0)
										{ break; }
									}
								}
								connector_start = connector_start_temp;
								
								if (match_connector)
								{
									
									// Initialize the AB ID
									String AB1_ID = "Unknown";
									String AB2_ID = "Unknown";
								
									// Initialize variables to store hamming distance for each AB1 and 2 barcode
									int[] temp1 = new int[ABarray.size()];
									int[] temp2 = new int[ABarray.size()];
									
									// Check if there is a frameshift, in order to locate AB2 correctly
									// Find the location of TAAAG in the found connector region
									int shift_j = line1.substring(connector_start, connector_start+connector.length()+2).indexOf("TAAAG"); // add 2 just in case there are 2 insertions
									if (shift_j == -1) // skip read if can't find TAAAG in the connector region
									{
										counter++;
										bad_connector_counter++;
//										bad_connector_counter_temp++;
										continue;
									}
									shift_j = 17 - shift_j; // number of bases to shift to the left
									
									// Found AB barcodes
									String AB1_found = line1.substring(connector_start-18, connector_start-18+8);
									String AB2_found = "";
									if ((connector_start+25-shift_j+8) <= line1.length()) // check if the read fully contains the AB2 ID
									{
										AB2_found = line1.substring(connector_start+25-shift_j, connector_start+25-shift_j+8);
									}
									else if ((connector_start+25-shift_j+8) == (line1.length()+1)) // the first 7 bases of the barcode is at the end of the read (in other words, some insertions)
									{
										AB2_found = line1.substring(connector_start+25-shift_j) + "N";
									}
									else // skip read if there are too many insertions
									{
										counter++;
										continue;
									}
									
									
									
									// Calculate Hamming distance only if there is at most 1 N in the found barcode
									if ((StringUtils.countMatches(AB1_found, "N") <= 1) && (StringUtils.countMatches(AB2_found, "N") <= 1))
									{
										// Counter for # of matches with 1 hamming distance
										int match_counter1 = 0, match_counter2 = 0;
										// Index for matches with 1 hamming distance
										int match_index1 = -1, match_index2 = -1;

										
										for (int i=0; (i<ABarray.size()) && (Objects.equals(AB1_ID,"Unknown") || Objects.equals(AB2_ID,"Unknown")) && ((match_counter1<=1) && (match_counter2<=1)); i++)
										{
											// Calculate Hamming distance
											temp1[i] = HammingDistanceCalculator(AB1_found, ABarray.get(i).get(1), true);
											temp2[i] = HammingDistanceCalculator(AB2_found, ABarray.get(i).get(1), true);
											
											// Allow early termination of for loop if found an exact match
											if (Objects.equals(AB1_ID,"Unknown"))
											{
												if (temp1[i] == 0)
												{
													AB1_ID = ABarray.get(i).get(0);
												}
												else if (temp1[i] == 1)
												{match_counter1++; match_index1 = i;}
											}
											
											if (Objects.equals(AB2_ID,"Unknown"))
											{
												if (temp2[i] == 0)
												{
													AB2_ID = ABarray.get(i).get(0);
												}
												else if (temp2[i] == 1)
												{match_counter2++; match_index2 = i;}
											}
											
										}
										
										// Find unambiguous match with 1 Hamming distance (ie, discard reads that have more than 1 matches with 1 hamming distance
										if (match_counter1 == 1) {AB1_ID = ABarray.get(match_index1).get(0);}
										if (match_counter2 == 1) {AB2_ID = ABarray.get(match_index2).get(0);}
										
										// There can be a deletion between the connector and the end of AB1 ID --> Check for this
										if ((match_counter1 == 0) && Objects.equals(AB1_ID,"Unknown"))
										{
											AB1_found = line1.substring(connector_start-18+1, connector_start-18+1+8);
											match_counter1 = 0;
											match_index1 = -1;
											
											for (int i=0; (i<ABarray.size()) && Objects.equals(AB1_ID,"Unknown") && (match_counter1<=1); i++)
											{
												// Calculate Hamming distance
												temp1[i] = HammingDistanceCalculator(AB1_found, ABarray.get(i).get(1), true);
												if (temp1[i] == 0)
												{
													AB1_ID = ABarray.get(i).get(0);
												}
												else if (temp1[i] == 1)
												{match_counter1++; match_index1 = i;}
											}
											
											if (match_counter1 == 1) {AB1_ID = ABarray.get(match_index1).get(0);}
										}
										
										if (!Objects.equals(AB1_ID,"Unknown") && !Objects.equals(AB2_ID,"Unknown"))
										{
											bwout.write(R1List_temp_array[1] + "," + line1.substring(1, 17) + "," + AB1_ID + "," + AB2_ID);
											bwout.newLine();
											PLA_counter++;
//											PLA_counter_temp++;
										}
										else
										{
											non_matching_AB_counter++;
//											non_matching_AB_counter_temp++;
										}
									}
									
									// Add the found AB barcodes to the Hash Multiset ABcounts
									ABcounts.add(AB1_found);
									ABcounts.add(AB2_found);
									
								}
								else
								{
									bad_connector_counter++;
//									bad_connector_counter_temp++;
								}
								
								
								
								if ((read_counter % 1_000_000) == 0)
								{
									System.out.printf("%s   ReadAlignmentSmartSeq   Processed %,15d reads   Elapsed time for last 1,000,000 reads: %ds%n",
											LocalDateTime.now().format(time_formatter), read_counter, (System.currentTimeMillis()-my_timer)/1000);
									my_timer = System.currentTimeMillis();
								}
							}
							
							
//							if ((((counter-1)/4)+1)>2000000) {System.out.printf("%s   %s   Processed %,d lines%n", LocalDateTime.now().format(time_formatter), counter); break;} // for testing purposes
						
							counter++;
						}
						
//						bwsum.write(R1List_temp_array[0]); bwsum.newLine();
//						bwsum.write("Number of reads with a valid PLA product: " + String.format("%,d", PLA_counter_temp)); bwsum.newLine();
//						bwsum.write("Number of reads discarded because of non-matching connector sequence: " + String.format("%,d",bad_connector_counter_temp)); bwsum.newLine();
//						bwsum.write("Number of reads discarded because of bad UMI: " + String.format("%,d",bad_UMI_counter_temp)); bwsum.newLine();
//						bwsum.write("Number of reads discarded because of non-matching antibody barcode: " + String.format("%,d",non_matching_AB_counter_temp)); bwsum.newLine();
//						bwsum.newLine();
						
					} catch (IOException e) {throw new IllegalArgumentException("Invalid file paths in R1List!");}
					
					
				}
				
				
				System.out.printf("%s   ReadAlignmentSmartSeq   Done: processed %,d reads%n", LocalDateTime.now().format(time_formatter), read_counter);
				System.out.printf("\tNumber of reads with a valid PLA product: %,15d%n", PLA_counter);
				System.out.println();
				
				// Write to summary file
				bwsum.write("ReadAlignmentSmartSeq: Finished at " + LocalDateTime.now().withNano(0) + ", processed " + String.format("%,d",read_counter) + " reads"); bwsum.newLine();
				bwsum.write("Total number of reads with a valid PLA product: " + String.format("%,d", PLA_counter)); bwsum.newLine();
				bwsum.write("Total number of reads discarded because of non-matching connector sequence: " + String.format("%,d",bad_connector_counter)); bwsum.newLine();
				bwsum.write("Total number of reads discarded because of bad UMI: " + String.format("%,d",bad_UMI_counter)); bwsum.newLine();
				bwsum.write("Total number of reads discarded because of non-matching antibody barcode: " + String.format("%,d",non_matching_AB_counter)); bwsum.newLine();

				// Add to the summary file the found AB barcodes
				bwsum.newLine();
				bwsum.write("Antibody barcode\tOccurrences"); bwsum.newLine();
				// Sort the HashMultiset by decreasing occurrences, and save the top 10 to an array
				String[] AB_sortedbycounts = Multisets.copyHighestCountFirst(ABcounts).elementSet().toArray(new String[0]);
				int temp_counter = 0;
				for (String i : AB_sortedbycounts)
				{
					// Break for loop once read count is lower than 100
					if (temp_counter >= 10)
					{break;}
					
					bwsum.write(i + String.format("\t%,10d", ABcounts.count(i)));
					bwsum.newLine();
					temp_counter++;
				}
				
			} catch (IOException e) {throw new IllegalArgumentException("Invalid file paths!");}
			break;
		}
		

		/**
		 * Cell barcode correction from aligned reads (ie, output of ReadAlignmentDropSeq)
		 * Use a list of reference cell barcodes to correct the cell barcodes of the reads
		 * 		I=... O=... SUMMARY=...
		 * 
		 * Input arguments:
		 * 		I: path to aligned reads (txt.gz format)
		 * 		O: path to store output file (txt.gz format)
		 * 		SUMMARY: directory to store summary files (txt format) (default is current working directory)
		 * 		CELL_BC_LIST: path to a comma-separated list of cell barcodes produced by drop-seq pipeline (rows are cell barcodes, column 0 is readcount, column 1 is the cell barcode sequence)
		 * 		^^^^^^^^^^^^ output of drop-seq tools' BAMTagHistogram function (txt.gz format)
		 * 		READCOUNT_CUTOFF: only keep the barcode sequence with at least this number of readcount (default is 1000)
		 * 		HEADER: whether the CELL_BC_LIST have header, which will be skipped (default is false)
		 * 
		 * Correction method: n-gram with Hamming distance <=1
		 * 		Split the query cell barcode into 2 halves at the first N
		 * 		Example: barcode "123N567N9" is split into "123" and "567N9"
		 * 		If there is no N, check for n-grams split in the middle: "123456" is split into "123" and "456"
		 * 		Then check for matching "123" and "5678" in the reference n-grams
		 * 		Look for perfectly matched cell barcodes, or cell barcodes that are unambiguously mismatched by 1 Hamming distance (ie, discard cell barcodes that have more than 1 matching cell barcodes with Hamming = 1) 
		 */	
		case "CellBarcodeCorrection":
		{
			
			// Parse the arguments
			String I = "", O = "", cell_BC_path = "";
			String SUMMARY = System.getProperty("user.dir") + File.separator + "CellBarcodeCorrection_summary.txt"; // default summary file directory
			int rc_cutoff = 1000;
			boolean skip_header = false;
			for (int i=1; i<args.length; i++)
			{
				String[] j = args[i].split("=");
				
				if (j.length < 2) {throw new java.lang.IllegalArgumentException("Whitespaces are not allowed between argument key specifier and the argument!");}
				
				switch (j[0])
				{
				case "I": I = j[1]; break;
				case "O": O = j[1]; break;
				case "SUMMARY": SUMMARY = j[1]; break;
				case "CELL_BC_LIST": cell_BC_path = j[1]; break;
				case "READCOUNT_CUTOFF": rc_cutoff = Integer.parseInt(j[1]); break;
				case "HEADER": skip_header = "true".equalsIgnoreCase(j[1]); break;
				default: throw new java.lang.IllegalArgumentException("Invalid argument key specifier!");
				}
			}
			
			try (
					BufferedReader brI = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(I)))); // aligned reads
					BufferedReader brBC = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(cell_BC_path)))); // cell barcode list
					BufferedWriter bwout = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(O)))); // output file
					BufferedWriter bwsum = new BufferedWriter(new FileWriter(SUMMARY)) // summary file
				)
			{
				
				// Write out the command line arguments
				for (String j : args)
				{
					bwsum.write(j + " "); bwsum.newLine();
					System.out.println(j);
				}
				bwsum.newLine();
				System.out.println();
				
				
				// Set up the counters for the summary files
				int PLA_counter = 0; // counter for the number of valid PLA reads
				int exact_counter = 0; // counter for the number of reads with exact match
				int ambiguous_counter = 0; // counter for the number of reads that matches ambiguously with reference cell barcodes (ie, more than 1 cell barcode with 1 hamming distance)
				int nomatch_counter = 0; // counter for the number of reads with no match
				int manyN_counter = 0; // counter for the number of reads with more than 1 N in the cell barcode region
				
				// Skip first line if required
				if (skip_header) {brBC.readLine();}
				
//				// Initialize a multimap of key-value pair: key is the sequence, and value is an ArrayList of the sequence's index in the reference barcode list
				ListMultimap<String, Integer> cellBCmap_1 = ArrayListMultimap.create(); // first segment of the barcode
				ListMultimap<String, Integer> cellBCmap_2 = ArrayListMultimap.create(); // second segment of the barcode
				
				// Read accepted cell barcode sequence into an ArrayList
				List<String> cellBCarray = new ArrayList<String>(); // ArrayList containing the cell barcode
				List<String> cellBCarray_1 = new ArrayList<String>(); // ArrayList containing the first half of the cell barcode (base 1-6)
				List<String> cellBCarray_2 = new ArrayList<String>(); // ArrayList containing the second half of the cell barcode (base 7-12)
				
				// Set up counters for reference cell barcodes
				int ref_counter = 0; // counter for the number of cell barcodes
				
				// Read the cell barcodes and build the multimap for the non-N case
				System.out.printf("%s   CellBarcodeCorrection   Start processing reference cell barcodes%n", LocalDateTime.now().format(time_formatter));
				String BCline;
				while ((BCline=brBC.readLine()) != null)
				{
					String[] values = BCline.split("\t");
					if (Integer.parseInt(values[0]) >= rc_cutoff)
					{
						cellBCarray.add(values[1]);
						// Split the reference barcode into 2 equal halves, and store them in the ArrayList
						cellBCarray_1.add(values[1].substring(0, values[1].length()/2));
						cellBCarray_2.add(values[1].substring(values[1].length()/2));
						
						// Store the 2 halves in the multimap
						cellBCmap_1.put(cellBCarray_1.get(cellBCarray_1.size()-1), ref_counter);
						cellBCmap_2.put(cellBCarray_2.get(cellBCarray_2.size()-1), ref_counter);
						ref_counter++;
					}
				}
				
				// Split the reference barcode at different positions for the N case
				// 		Example: barcode "abcdef", N can be at any position
				// Initialize a list of multimaps of key-value pair
				// Each multimap is the n-gram and its index in the original reference list, split at a specific location
				//		Example: the 2nd multimap in the list contains the n-gram split at the 2nd base: "aNcdef"
				// If N is at the beginning or the end of the cell barcode, ignore N, and split at the middle: "Nbcdef" is split into "bc" and "def", "abcdeN" is split into "abc" and "de"
				List<ListMultimap<String, Integer>> cellBCmaplist_1 = new ArrayList<>();
				List<ListMultimap<String, Integer>> cellBCmaplist_2 = new ArrayList<>();
				
				// Initialize an ArrayList of ArrayList
				// Each sub-ArrayList is all the segments of the reference cell barcode (non-collapsed), split at a specific location
				List<ArrayList<String>> cellBCarraylist_1 = new ArrayList<>();
				List<ArrayList<String>> cellBCarraylist_2 = new ArrayList<>();
				
				System.out.printf("%s   CellBarcodeCorrection   Start building maps of reference cell barcodes%n", LocalDateTime.now().format(time_formatter));
				int barcode_length = cellBCarray.get(0).length();
				for (int i=0; i<barcode_length; i++) // loop through each possible splitting position
				{
					// Temporary variables
					ListMultimap<String, Integer> tempmap_1 = ArrayListMultimap.create();
					ListMultimap<String, Integer> tempmap_2 = ArrayListMultimap.create();
					ArrayList<String> temparray_1 = new ArrayList<>();
					ArrayList<String> temparray_2 = new ArrayList<>();
					
					// Specify start index (inclusive) of the first n-gram, and end index (exclusive) of the second n-gram
					int temp_start = 0;
					int temp_end = barcode_length;
					// Specify split position
					int temp_split_1 = i;
					int temp_split_2 = i + 1;
					if (i == 0)
					{
						temp_start = 1;
						temp_split_1 = barcode_length/2;
						temp_split_2 = barcode_length/2;
					}
					else if (i == (barcode_length-1))
					{
						temp_end = barcode_length - 1;
						temp_split_1 = barcode_length/2;
						temp_split_2 = barcode_length/2;
					}
					
					for (int j=0; j<cellBCarray.size(); j++) // loop through reference cell barcodes
					{
						String temp_1 = cellBCarray.get(j).substring(temp_start, temp_split_1);
						String temp_2 = cellBCarray.get(j).substring(temp_split_2, temp_end);
						tempmap_1.put(temp_1, j);
						tempmap_2.put(temp_2, j);
						
						temparray_1.add(temp_1);
						temparray_2.add(temp_2);
					}
					
					cellBCmaplist_1.add(tempmap_1);
					cellBCmaplist_2.add(tempmap_2);
					cellBCarraylist_1.add(temparray_1);
					cellBCarraylist_2.add(temparray_2);
				}
				
				System.out.printf("%s   CellBarcodeCorrection   Finish reading %,d cell barcodes%n", LocalDateTime.now().format(time_formatter), ref_counter);
				
				
				// n-gram method: takes into account N
				String line;
				int counter = 0;
				long my_timer = System.currentTimeMillis();
				System.out.printf("%s   Start cell barcode correction%n", LocalDateTime.now().format(time_formatter));
				bwsum.write("Start cell barcode correction at " + LocalDateTime.now()); bwsum.newLine();
				while ((line=brI.readLine()) != null)
				{
						
					// Time stamp (put here because there are continue statements below)
					if ((counter>0) && (counter % 1000000 == 0))
					{
						System.out.printf("%s   CellBarcodeCorrection   Processed %,15d records   Elapsed time for last 1,000,000 records: %ds%n",
								LocalDateTime.now().format(time_formatter), counter, (System.currentTimeMillis()-my_timer)/1000);
						my_timer = System.currentTimeMillis();
					}
					
					// Read the cell barcodes, UMI and AB IDs
					List<String> values = new ArrayList<String>(Arrays.asList(line.split(",")));
					
					// counter for how many barcodes have hamming distance <= 1, discard cell barcode if counter >=2
					int match_counter = 0;
					
					// Case 1: query barcode doesn't contain N
					if (!(values.get(0).contains("N")))
					{
						// Split cell barcodes into half, and compare each to the reference cell barcodes
						String testBC_1 = values.get(0).substring(0, values.get(0).length()/2);
						String testBC_2 = values.get(0).substring(values.get(0).length()/2);
						
						// Check if the first half has a match
						if (cellBCmap_1.containsKey(testBC_1))
						{
							// Look for matching second half
							for (int i : cellBCmap_1.get(testBC_1))
							{
								if (Objects.equals(testBC_2, cellBCarray_2.get(i)))
								{
									exact_counter++;
									match_counter = 1;
									break;
								}
								else if (HammingDistanceCalculator(testBC_2,cellBCarray_2.get(i)) == 1)
								{
									match_counter++;
									values.set(0, cellBCarray.get(i));
								}
							}
						}
						// Check if the second half has a match
						else if (cellBCmap_2.containsKey(testBC_2))
						{
							// Look for matching first half with 1 hamming distance
							for (int i : cellBCmap_2.get(testBC_2))
							{
								if (HammingDistanceCalculator(testBC_1, cellBCarray_1.get(i)) == 1)
								{
									match_counter++;
									values.set(0, cellBCarray_1.get(i) + testBC_2);
								}
								if (match_counter >= 2)
								{
									break;
								}
								
							}
						}
					}
					
					// Case 2: query barcode does contain N
					else
					{
						// Allow up to 1 N's only
						if (StringUtils.countMatches(values.get(0), "N") > 1)
						{
							manyN_counter++;
							counter++;
							continue;
						}
						
						int split_position = values.get(0).indexOf("N"); // position of N
						
						// Split query cell barcode into n-gram
						String testBC_1, testBC_2;
						if (split_position == 0)
						{
							testBC_1 = values.get(0).substring(1, values.get(0).length()/2);
							testBC_2 = values.get(0).substring(values.get(0).length()/2 );
						}
						else if (split_position == (values.get(0).length()-1))
						{
							testBC_1 = values.get(0).substring(0, values.get(0).length()/2);
							testBC_2 = values.get(0).substring(values.get(0).length()/2, values.get(0).length() - 1);
						}
						else
						{
							testBC_1 = values.get(0).substring(0, split_position);
							testBC_2 = values.get(0).substring(split_position + 1);
						}
						
						
						// Check if the first half has a match
						if (cellBCmaplist_1.get(split_position).containsKey(testBC_1))
						{
							for (int i : cellBCmaplist_1.get(split_position).get(testBC_1))
							{
								if (Objects.equals(testBC_2, cellBCarraylist_2.get(split_position).get(i)))
								{
									exact_counter++;
									match_counter = 1;
									values.set(0, cellBCarray.get(i));
									break;
								}
								else if (HammingDistanceCalculator(testBC_2, cellBCarraylist_2.get(split_position).get(i)) == 1)
								{
									match_counter++;
									values.set(0, cellBCarray.get(i));
								}
							}
						}
						
						else if (cellBCmaplist_2.get(split_position).containsKey(testBC_2))
						{
							for (int i : cellBCmaplist_2.get(split_position).get(testBC_2))
							{
								if (HammingDistanceCalculator(testBC_1, cellBCarraylist_1.get(split_position).get(i)) == 1)
								{
									match_counter++;
									values.set(0, cellBCarray.get(i));
								}
								
								if (match_counter >= 2)
								{
									break;
								}
							}
						}
					}
					
								
					if (match_counter == 1)
					{
						PLA_counter++;
						bwout.write(String.join(",", values));
						bwout.newLine();
					}
					else if (match_counter == 0) // no matching barcode is found
					{
						nomatch_counter++;
					}
					else if (match_counter >= 2)
					{
						ambiguous_counter++;
					}
					
					counter++;
					
				}
				
				System.out.printf("%s   CellBarcodeCorrection   Done%n", LocalDateTime.now().format(time_formatter));
				System.out.printf("\tNumber of valid PLA products: %,15d%n", PLA_counter);
				
				// Write to summary file
				bwsum.write("CellBarcodeCorrection: Finished at " + LocalDateTime.now().withNano(0) + ", processed " + String.format("%,d",counter) + " records"); bwsum.newLine();
				bwsum.write("Number of reference cell barcodes " + String.format("%,d", cellBCarray.size())); bwsum.newLine();
				bwsum.write("Number of accepted PLA products: " + String.format("%,d", PLA_counter)); bwsum.newLine();
				bwsum.write("Number of exact matches: " + String.format("%,d",exact_counter)); bwsum.newLine();
				bwsum.write("Number of discarded ambiguous reads: " + String.format("%,d",ambiguous_counter)); bwsum.newLine();
				bwsum.write("Number of discarded non-matching reads: " + String.format("%,d",nomatch_counter)); bwsum.newLine();
				bwsum.write("Number of discarded cell barcodes with more than 1 N base: " + String.format("%,d", manyN_counter)); bwsum.newLine();
				
			} catch (IOException e) {throw new IllegalArgumentException("Invalid file paths!");}
			break;
		}

		
		/**
		 * Perform cell barcode correction from aligned reads (ie, output of ReadAlignmentDropSeq)
		 * Use a list of reference cell barcodes to correct the cell barcodes of the reads
		 * Does NOT use any reference cell barcodes (while CellBarcodeCorrection does)
		 * 		I=... O=... SUMMARY=...
		 * 
		 * Input arguments:
		 * 		I: path to aligned/cell barcode-corrected file (txt.gz format)
		 * 		O: path to output UMI merged file (txt.gz format)
		 * 		SUMMARY: path to store summary file (txt format)
		 * 
		 * All detected cell barcodes are merged as following:
		 * 		Exact matches
		 * 		Match with 1 hamming distance to another UMI that has a higher read count
		 * 	
		 * The merged cell barcodes are then used as the reference cell barcode list, and proceed as CellBarcodeCorrection
		 */
		
		case "CellBarcodeMerging":
		{
			
			// Parse the arguments
			String I = "", O = "";
			String SUMMARY = System.getProperty("user.dir") + File.separator + "CellBarcodeCorrection_summary.txt"; // default summary file directory
			for (int i=1; i<args.length; i++)
			{
				String[] j = args[i].split("=");
				
				if (j.length < 2) {throw new java.lang.IllegalArgumentException("Whitespaces are not allowed between argument key specifier and the argument!");}
				
				switch (j[0])
				{
				case "I": I = j[1]; break;
				case "O": O = j[1]; break;
				case "SUMMARY": SUMMARY = j[1]; break;
				default: throw new java.lang.IllegalArgumentException("Invalid argument key specifier!");
				}
			}
			
			// Write out the command line arguments
			for (String j : args)
			{
				System.out.println(j);
			}
			System.out.println();
			
			// Create ArrayList of reference cell barcodes, obtained from pre-processing the detected cell barcodes
			List<String> refBC = new ArrayList<String>();
			
			// Read the file once to collect the cell barcodes
			try (BufferedReader brI = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(I)))))
			{
				// Create Hash Multiset to get counts of each cell barcode
				Multiset<String> BC_readcounts = HashMultiset.create();
				
				// Read cell barcode sequence into an ArrayList
				// Ignore cell barcodes containing N
				System.out.printf("%s   CellBarcodeMerging   Collecting detected cell barcodes%n", LocalDateTime.now().format(time_formatter));
				String BCline;
				while ((BCline = brI.readLine()) != null)
				{
					String[] values = BCline.split(",");
					if (!(values[0].contains("N")))
					{
						BC_readcounts.add(values[0]);
					}
				}
				
				// Sort the HashMultiset by decreasing read counts, and save to an array
				String[] BC_sortedbycounts = Multisets.copyHighestCountFirst(BC_readcounts).elementSet().toArray(new String[0]);
				
				// Split accepted cell barcode sequence into two halves, and store in two ArrayList's
				List<String> cellBCarray_1 = new ArrayList<String>(); // ArrayList containing the first half of the cell barcode (base 1-6)
				List<String> cellBCarray_2 = new ArrayList<String>(); // ArrayList containing the second half of the cell barcode (base 7-12)
				
				// Initialize a multimap of key-value pair: key is the sequence, and value is an ArrayList of the sequence's indices in BC_sortedbycounts
				// The index in the ArrayList value will be used to index cellBCarray_1 and cellBCarray_2
				ListMultimap<String, Integer> cellBCmap_1 = ArrayListMultimap.create(); // first half of the barcode
				ListMultimap<String, Integer> cellBCmap_2 = ArrayListMultimap.create(); // second half of the barcode
				
				System.out.printf("%s   CellBarcodeMerging   Merginng collected cell barcodes%n", LocalDateTime.now().format(time_formatter));
				
				// Merge cell barcodes: convert a cell barcode to another that is at 1 hamming distance and has a higher read count
				// Split detected cell barcodes into two halves, and use n-gram to merge cell barcodes
				int BC_index = 0;
				int BC_length = BC_sortedbycounts[0].length();
				for (String BC_i : BC_sortedbycounts)
				{					
					// Update cellBCarray_1 and _2
					cellBCarray_1.add(BC_i.substring(0, BC_length/2));
					cellBCarray_2.add(BC_i.substring(BC_length/2));
					
					cellBCmap_1.put(cellBCarray_1.get(cellBCarray_1.size()-1), BC_index);
					cellBCmap_2.put(cellBCarray_2.get(cellBCarray_2.size()-1), BC_index);
					BC_index++;
				}
				
				// LinkedHashSet to store unique cell barcodes, from which converted barcodes have been removed
				Set<String> BC_updated = new LinkedHashSet<>(Arrays.asList(BC_sortedbycounts)); // use LinkedHashSet to keep barcode ordered by read count
				// HashSet to store barcodes that have been removed
				Set<String> BC_removed = new HashSet<>();
				
				// Iterate from lowest read count barcode to highest, and remove the one that is 1 hamming distance away from another barcode with the highest read count
				for (int j=(BC_sortedbycounts.length-1); j>=0; j--)
				{
					
					if (BC_removed.contains(BC_sortedbycounts[j])) {break;}
					
					// Get the first and second half of the barcode
					String temp_1 = BC_sortedbycounts[j].substring(0, BC_length/2);
					String temp_2 = BC_sortedbycounts[j].substring(BC_length/2);
					
					// Get all matches that are 1 hamming distance away
					// Check if the first half has a match
					if (cellBCmap_1.containsKey(temp_1))
					{
						// Look for matching second half
						for (int i : cellBCmap_1.get(temp_1))
						{
							if (HammingDistanceCalculator(temp_2, cellBCarray_2.get(i)) == 1)
							{
								// A check to prevent "correct" cell barcodes from being removed
								if (BC_removed.contains(BC_sortedbycounts[i])) {break;}
								
								// Remove match barcode
								BC_updated.remove(BC_sortedbycounts[j]);
								BC_removed.add(BC_sortedbycounts[j]);
								break;
							}
						}
					}
					// Check if the second half has a match
					else if (cellBCmap_2.containsKey(temp_2))
					{
						// Look for matching first half with 1 hamming distance
						for (int i : cellBCmap_2.get(temp_2))
						{
							if (HammingDistanceCalculator(temp_1, cellBCarray_1.get(i)) == 1)
							{
								// A check to prevent "correct" cell barcodes from being removed
								if (BC_removed.contains(BC_sortedbycounts[i])) {break;}
								
								// Remove match barcode
								BC_updated.remove(BC_sortedbycounts[j]);
								BC_removed.add(BC_sortedbycounts[j]);
								break;
							}		
						}
					}
				}
				
				// Export the reference cell barcodes
				refBC.addAll(0, BC_updated);
				
				System.out.printf("%s   CellBarcodeMerging   Done merginng collected cell barcodes%n", LocalDateTime.now().format(time_formatter));
				
			} catch (IOException e) {throw new IllegalArgumentException("Invalid input file path!");}
			
			try (
					BufferedReader brI = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(I)))); // aligned reads
					BufferedWriter bwout = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(O)))); // output file
					BufferedWriter bwsum = new BufferedWriter(new FileWriter(SUMMARY)) // summary file
				)
			{
				
				// Write out the command line arguments
				for (String j : args)
				{
					bwsum.write(j + " "); bwsum.newLine();
				}
				bwsum.newLine();
				
				
				// Set up the counters for the summary files
				int PLA_counter = 0; // counter for the number of valid PLA reads
				int exact_counter = 0; // counter for the number of reads with exact match
				int ambiguous_counter = 0; // counter for the number of reads that matches ambiguously with reference cell barcodes (ie, more than 1 cell barcode with 1 hamming distance)
				int nomatch_counter = 0; // counter for the number of reads with no match
				int manyN_counter = 0; // counter for the number of reads with more than 1 N in the cell barcode region
				
				System.out.printf("%s   CellBarcodeMerging   Start processing reference cell barcodes%n", LocalDateTime.now().format(time_formatter));
				
				// Read accepted cell barcode sequence into an ArrayList
				List<String> cellBCarray = new ArrayList<String>(); // ArrayList containing the cell barcode
				List<String> cellBCarray_1 = new ArrayList<String>(); // ArrayList containing the first half of the cell barcode (base 1-6)
				List<String> cellBCarray_2 = new ArrayList<String>(); // ArrayList containing the second half of the cell barcode (base 7-12)
				
				// Initialize a multimap of key-value pair: key is the sequence, and value is an ArrayList of the sequence's index in the reference barcode list
				// The index in the ArrayList value will be used to index cellBCarray_1 and cellBCarray_2
				ListMultimap<String, Integer> cellBCmap_1 = ArrayListMultimap.create(); // first half of the barcode
				ListMultimap<String, Integer> cellBCmap_2 = ArrayListMultimap.create(); // second half of the barcode
				
				// Set up counters for reference cell barcodes
				int ref_counter = 0; // counter for the number of cell barcodes
				
				// Read the cell barcodes and build the multimap for the non-N case
				int barcode_length = refBC.get(0).length();
				for (String BC_i : refBC)
				{
					cellBCarray.add(BC_i);
					// Split the reference barcode into 2 equal halves, and store them in the ArrayList
					cellBCarray_1.add(BC_i.substring(0, barcode_length/2));
					cellBCarray_2.add(BC_i.substring(barcode_length/2));
					
					// Store the 2 halves in the multimap
					cellBCmap_1.put(cellBCarray_1.get(cellBCarray_1.size()-1), ref_counter);
					cellBCmap_2.put(cellBCarray_2.get(cellBCarray_2.size()-1), ref_counter);
					ref_counter++;
				}
				
				// Split each reference barcode at different positions for the N case
				// 		Example: barcode "abcdef", N can be at any position
				// Initialize a list of multimaps of key-value pair
				// Each multimap is the n-gram and its index in the original reference list, split at a specific location
				//		Example: the 2nd multimap in the list contains the n-gram split at the 2nd base: "aNcdef"
				// If N is at the beginning or the end of the cell barcode, ignore N, and split at the middle: "Nbcdef" is split into "bc" and "def", "abcdeN" is split into "abc" and "de"
				List<ListMultimap<String, Integer>> cellBCmaplist_1 = new ArrayList<>();
				List<ListMultimap<String, Integer>> cellBCmaplist_2 = new ArrayList<>();
				
				// Initialize an ArrayList of ArrayList
				// Each sub-ArrayList is all the segments of the reference cell barcode (non-collapsed), split at a specific location
				List<ArrayList<String>> cellBCarraylist_1 = new ArrayList<>();
				List<ArrayList<String>> cellBCarraylist_2 = new ArrayList<>();
				
				System.out.printf("%s   CellBarcodeMerging   Start building maps of reference cell barcodes%n", LocalDateTime.now().format(time_formatter));
				
				for (int i=0; i<barcode_length; i++) // loop through each possible splitting position
				{
					// Temporary variables
					ListMultimap<String, Integer> tempmap_1 = ArrayListMultimap.create();
					ListMultimap<String, Integer> tempmap_2 = ArrayListMultimap.create();
					ArrayList<String> temparray_1 = new ArrayList<>();
					ArrayList<String> temparray_2 = new ArrayList<>();
					
					// Specify start index (inclusive) of the first n-gram, and end index (exclusive) of the second n-gram
					int temp_start = 0;
					int temp_end = barcode_length;
					// Specify split position
					int temp_split_1 = i;
					int temp_split_2 = i + 1;
					if (i == 0)
					{
						temp_start = 1;
						temp_split_1 = barcode_length/2;
						temp_split_2 = barcode_length/2;
					}
					else if (i == (barcode_length-1))
					{
						temp_end = barcode_length - 1;
						temp_split_1 = barcode_length/2;
						temp_split_2 = barcode_length/2;
					}
					
					for (int j=0; j<cellBCarray.size(); j++) // loop through reference cell barcodes
					{
						String temp_1 = cellBCarray.get(j).substring(temp_start, temp_split_1);
						String temp_2 = cellBCarray.get(j).substring(temp_split_2, temp_end);
						tempmap_1.put(temp_1, j);
						tempmap_2.put(temp_2, j);
						
						temparray_1.add(temp_1);
						temparray_2.add(temp_2);
					}
					
					cellBCmaplist_1.add(tempmap_1);
					cellBCmaplist_2.add(tempmap_2);
					cellBCarraylist_1.add(temparray_1);
					cellBCarraylist_2.add(temparray_2);
				}
				
				
				// n-gram method: takes into account N
				String line;
				int counter = 0;
				long my_timer = System.currentTimeMillis();
				System.out.printf("%s   CellBarcodeMerging   Start cell barcode correction%n", LocalDateTime.now().format(time_formatter));
				bwsum.write("Start cell barcode correction at " + LocalDateTime.now()); bwsum.newLine();
				while ((line=brI.readLine()) != null)
				{
						
					// Time stamp (put here because there are continue statements below)
					if ((counter>0) && (counter % 1000000 == 0))
					{
						System.out.printf("%s   CellBarcodeMerging   Processed %,15d records   Elapsed time for last 1,000,000 records: %ds%n",
								LocalDateTime.now().format(time_formatter), counter, (System.currentTimeMillis()-my_timer)/1000);
						my_timer = System.currentTimeMillis();
					}
					
					// Read the cell barcodes, UMI and AB IDs
					List<String> values = new ArrayList<String>(Arrays.asList(line.split(",")));
					
					// counter for how many barcodes have hamming distance <= 1, discard cell barcode if counter >=2
					int match_counter = 0;
					
					// Case 1: query barcode doesn't contain N
					if (!(values.get(0).contains("N")))
					{
						// Split cell barcodes into half, and compare each to the reference cell barcodes
						String testBC_1 = values.get(0).substring(0, barcode_length/2);
						String testBC_2 = values.get(0).substring(barcode_length/2);
						
						// Check if the first half has a match
						if (cellBCmap_1.containsKey(testBC_1))
						{
							// Look for matching second half
							for (int i : cellBCmap_1.get(testBC_1))
							{
								if (Objects.equals(testBC_2, cellBCarray_2.get(i)))
								{
									exact_counter++;
									match_counter = 1;
									break;
								}
								else if (HammingDistanceCalculator(testBC_2,cellBCarray_2.get(i)) == 1)
								{
									match_counter++;
									values.set(0, cellBCarray.get(i));
								}
							}
						}
						// Check if the second half has a match
						else if (cellBCmap_2.containsKey(testBC_2))
						{
							// Look for matching first half with 1 hamming distance
							for (int i : cellBCmap_2.get(testBC_2))
							{
								if (HammingDistanceCalculator(testBC_1, cellBCarray_1.get(i)) == 1)
								{
									match_counter++;
									values.set(0, cellBCarray_1.get(i) + testBC_2);
								}
								if (match_counter >= 2)
								{
									break;
								}
								
							}
						}
					}
					
					// Case 2: query barcode does contain N
					else
					{
						// Allow up to 1 N's only
						if (StringUtils.countMatches(values.get(0), "N") > 1)
						{
							manyN_counter++;
							counter++;
							continue;
						}
						
						int split_position = values.get(0).indexOf("N"); // position of N
						
						// Split query cell barcode into n-gram
						String testBC_1, testBC_2;
						if (split_position == 0)
						{
							testBC_1 = values.get(0).substring(1, barcode_length/2);
							testBC_2 = values.get(0).substring(barcode_length/2);
						}
						else if (split_position == (barcode_length-1))
						{
							testBC_1 = values.get(0).substring(0, barcode_length/2);
							testBC_2 = values.get(0).substring(barcode_length/2, barcode_length - 1);
						}
						else
						{
							testBC_1 = values.get(0).substring(0, split_position);
							testBC_2 = values.get(0).substring(split_position + 1);
						}
						
						
						// Check if the first half has a match
						if (cellBCmaplist_1.get(split_position).containsKey(testBC_1))
						{
							for (int i : cellBCmaplist_1.get(split_position).get(testBC_1))
							{
								if (Objects.equals(testBC_2, cellBCarraylist_2.get(split_position).get(i)))
								{
									exact_counter++;
									match_counter = 1;
									values.set(0, cellBCarray.get(i));
									break;
								}
								else if (HammingDistanceCalculator(testBC_2, cellBCarraylist_2.get(split_position).get(i)) == 1)
								{
									match_counter++;
									values.set(0, cellBCarray.get(i));
								}
							}
						}
						
						else if (cellBCmaplist_2.get(split_position).containsKey(testBC_2))
						{
							for (int i : cellBCmaplist_2.get(split_position).get(testBC_2))
							{
								if (HammingDistanceCalculator(testBC_1, cellBCarraylist_1.get(split_position).get(i)) == 1)
								{
									match_counter++;
									values.set(0, cellBCarray.get(i));
								}
								
								if (match_counter >= 2)
								{
									break;
								}
							}
						}
					}
					
								
					if (match_counter == 1)
					{
						PLA_counter++;
						bwout.write(String.join(",", values));
						bwout.newLine();
					}
					else if (match_counter == 0) // no matching barcode is found
					{
						nomatch_counter++;
					}
					else if (match_counter >= 2)
					{
						ambiguous_counter++;
					}
					
					counter++;
					
				}
				
				System.out.printf("%s   CellBarcodeMerging   Done%n", LocalDateTime.now().format(time_formatter));
				System.out.printf("\tNumber of valid PLA products: %,15d%n", PLA_counter);
				
				// Write to summary file
				bwsum.write("CellBarcodeMerging: Finished at " + LocalDateTime.now().withNano(0) + ", processed " + String.format("%,d",counter) + " records"); bwsum.newLine();
				bwsum.write("Number of reference cell barcodes " + String.format("%,d",refBC.size())); bwsum.newLine();
				bwsum.write("Number of accepted PLA products: " + String.format("%,d", PLA_counter)); bwsum.newLine();
				bwsum.write("Number of exact matches: " + String.format("%,d",exact_counter)); bwsum.newLine();
				bwsum.write("Number of discarded ambiguous reads: " + String.format("%,d",ambiguous_counter)); bwsum.newLine();
				bwsum.write("Number of discarded non-matching reads: " + String.format("%,d",nomatch_counter)); bwsum.newLine();
				bwsum.write("Number of discarded cell barcodes with more than 1 N base: " + String.format("%,d", manyN_counter)); bwsum.newLine();
				
			} catch (IOException e) {throw new IllegalArgumentException("Invalid file paths!");}
			break;
		}
		
		
		/**
		 * Check if any UMIs map to more than 1 PLA product
		 * 		I=... O=...
		 * Input arguments:
		 * 		I: path to aligned/cell barcode-corrected file (txt.gz format)
		 * 		O: path to output multi-mapping UMIs (tab-separated, txt.gz format)
		 * 
		 * Output:
		 * 		Each row contains: cell barcode_UMI [PLA product 1, PLA product 2, ...]
		 * 		PLA product 1, 2, ... all have the same cell barcode and UMI
		 */
		
		case "CheckUMIMapping":
		{
			// Parse the arguments
			String I = "", O = "";
			for (int i=1; i<args.length; i++)
			{
				String[] j = args[i].split("=");
				
				if (j.length < 2) {throw new java.lang.IllegalArgumentException("Whitespaces are not allowed between argument key specifier and the argument!");}
				
				switch (j[0])
				{
				case "I": I = j[1]; break;
				case "O": O = j[1]; break;
				default: throw new java.lang.IllegalArgumentException("Invalid argument key specifier!");
				}
			}
			
			try (
					BufferedReader brI = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(I)))); // aligned reads
					BufferedWriter bwout = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(O)))) // output file
				  )
			{
				// Write out the command line arguments
				for (String j : args)
				{
					System.out.println(j);
				}
				System.out.printf("%n");
				
				// Initialize a SetMultimap
				// keys: unique cell barcode-UMI combinations (sorted alphabetically)
				// values: set of corresponding unique PLA products
				SetMultimap<String, String> UMI_multimap = MultimapBuilder.treeKeys().hashSetValues().build();
				
				// Time stamp
				long my_timer = System.currentTimeMillis();
				System.out.printf("%s   CheckUMIMapping   Start%n", LocalDateTime.now().format(time_formatter));
				
				// Add the aligned reads to the multimap
				String line;
				while ((line = brI.readLine()) != null)
				{
					String[] values = line.split(",");
					UMI_multimap.put(values[0] + "_" + values[1], values[2] + ":" + values[3]);
				}
				
				// Export duplicated UMI mappings
				int counter = 0;
				Set<String> temp_values = new HashSet<>();
				for (String i : UMI_multimap.keySet())
				{
					temp_values = UMI_multimap.get(i);
					if (temp_values.size() > 1)
					{
						bwout.write(i + "\t" + temp_values);
						bwout.newLine();
					}
					counter++;
					
					if ((counter % 100_000) == 0)
					{
						System.out.printf("%s   CheckUMIMapping   Processed %,15d UMIs   Elapsed time for last 100,000 UMIs: %ds%n",
								LocalDateTime.now().format(time_formatter), (counter-1)/4+1, (System.currentTimeMillis()-my_timer)/1000);
						my_timer = System.currentTimeMillis();
					}
				}
				
				// Time stamp
				System.out.printf("%s   CheckUMIMapping   Done: processed %,d UMIs%n", LocalDateTime.now().format(time_formatter), counter);
				
			} catch (IOException e) { throw new IllegalArgumentException("Invalid file paths!");}
			
			break;
		}
		
		
		/**
		 * Perform UMI merging on cell barcode-corrected PLA products
		 * 		I=... O=... SUMMARY=...
		 * Input arguments:
		 * 		I: path to aligned/cell barcode-corrected file (txt.gz format)
		 * 		O: path to output UMI merged file (txt.gz format)
		 * 		SUMMARY: path to store summary file
		 * 
		 * Merge the following cases:
		 * 		Exact matches
		 * 		Match with 1 hamming distance to another UMI that has a higher read count (or the highest read count if there're more than 1 matches)
		 */
		
		case "UMIMerging":
		{
			// Parse the arguments
			String I = "", O = "";
			String SUMMARY = System.getProperty("user.dir") + File.separator + "UMIMerging_summary.txt"; // default summary file directory
			for (int i=1; i<args.length; i++)
			{
				String[] j = args[i].split("=");
				
				if (j.length < 2) {throw new java.lang.IllegalArgumentException("Whitespaces are not allowed between argument key specifier and the argument!");}
				
				switch (j[0])
				{
				case "I": I = j[1]; break;
				case "O": O = j[1]; break;
				case "SUMMARY": SUMMARY = j[1]; break;
				default: throw new java.lang.IllegalArgumentException("Invalid argument key specifier!");
				}
			}
			
			try (
					BufferedReader brI = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(I)))); // aligned reads
					BufferedWriter bwout = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(O)))); // output file
					BufferedWriter bwsum = new BufferedWriter(new FileWriter(SUMMARY)) // summary file
				  )
			{

				// Write out the command line arguments
				for (String j : args)
				{
					bwsum.write(j); bwsum.newLine();
					System.out.println(j);
				}
				bwsum.newLine();
				System.out.printf("%n");
				
				// Each PLA product is described by 2 features: [cell barcode/AB1 ID/AB2 ID], and [UMI]. We shall call these feature 1 and 2, respectively
				// Initialize ArrayListMultimap: key is unique feature 1's, value is its corresponding feature 2's
				ListMultimap<String, String> PLAproductmap = ArrayListMultimap.create();
				
				// Read the file
				System.out.printf("%s   UMIMerging   Start pre-processing records%n", LocalDateTime.now().format(time_formatter));
				bwsum.write("Start UMIMerging at " + LocalDateTime.now().withNano(0)); bwsum.newLine();
				long my_timer = System.currentTimeMillis();
				String line;
				int record_counter = 0; // for time keeping
				while ((line = brI.readLine()) != null)
				{
					String[] values = line.split(",");
					PLAproductmap.put(values[0] + "," + values[2] + "," + values[3], values[1]); // add feature 1 and feature 2
					record_counter++;
					if ((record_counter % 1_000_000) == 0)
					{
						System.out.printf("%s   UMIMerging   Processed %,15d records   Elapsed time for last 1,000,000 records: %ds%n",
								LocalDateTime.now().format(time_formatter), record_counter, (System.currentTimeMillis()-my_timer)/1000);
						my_timer = System.currentTimeMillis();
					}
				}
				
				System.out.printf("%s   UMIMerging   Finished pre-processing%n", LocalDateTime.now().format(time_formatter));
				System.out.printf("\tNumber of unique cell barcode-PLA pair combinations: %,15d%n", PLAproductmap.keySet().size());
				
				// Loop through each unique feature 1, and merge UMI from the same feature 1 that are within 1 hamming distance
				int counter = 0; // total record counter
				int unique_counter = 0; // unique PLA product counter
				my_timer = System.currentTimeMillis();
				System.out.printf("%s   UMIMerging   Start merging UMI%n", LocalDateTime.now().format(time_formatter));
				for (String str_i : PLAproductmap.keySet())
				{

					// Hash Multiset to store unique UMIs and their number of occurrences
					Multiset<String> temp_UMIcounts = HashMultiset.create();
					
					// Get the UMIs and their read counts associated with this feature 1
					for (String j : PLAproductmap.get(str_i))
					{
						temp_UMIcounts.add(j);
					}
					
					// Hashset to store unique UMIs, from which duplicated UMIs have been removed
					Set<String> temp_UMI_updated = new HashSet<>(temp_UMIcounts.elementSet());
					
					// Sort the HashMultiset by decreasing occurrences, and save to an array
					String[] temp_UMI_sortedbycounts = Multisets.copyHighestCountFirst(temp_UMIcounts).elementSet().toArray(new String[0]);
					
					// Iterate through the unique UMIs in order of increasing occurrences, and convert UMIs to UMIs that are within 1 hamming distance and have higher read counts
					// ****Slow step**** probably faster with BK-Tree //
					for (int j=(temp_UMI_sortedbycounts.length-1); j>=0; j--)
					{
						// Temporary set for iteration
						Set<String> temp_UMI = new HashSet<>(temp_UMI_updated);
						for (String str_j : temp_UMI)
						{
							if (HammingDistanceCalculator(str_j, temp_UMI_sortedbycounts[j], true) == 1)
							{
								temp_UMI_updated.remove(temp_UMI_sortedbycounts[j]);
								break;
							}
						}
					}
					
					// Export the unique UMI reads
					String[] values = str_i.split(",");
					for (String str_j : temp_UMI_updated)
					{
						bwout.write(values[0] + "," + str_j + "," + values[1] + "," + values[2]); bwout.newLine();
						unique_counter++;
					}
					
					// Time stamp
					counter++;
					if ((counter % 1_000) == 0)
					{
						System.out.printf("%s   UMIMerging   Processed %,12d records   Elapsed time for last 1,000 records: %ds%n",
								LocalDateTime.now().format(time_formatter), counter, (System.currentTimeMillis()-my_timer)/1000);
						my_timer = System.currentTimeMillis();
					}
				}
				
				System.out.printf("%s   UMIMerging   Done%n", LocalDateTime.now().format(time_formatter));
				System.out.printf("\tNumber of unique PLA products: %,15d%n", unique_counter);
				
				// Write to summary file
				bwsum.write("UMIMerging: Finished at " + LocalDateTime.now().withNano(0) + ", processed " + String.format("%,d",counter) + " unique cell barcode-PLA pair combinations"); bwsum.newLine();
				bwsum.write("Number of unique PLA products: " + String.format("%,d", unique_counter)); bwsum.newLine(); bwsum.newLine();

				
				
			} catch (IOException e) { throw new IllegalArgumentException("Invalid file paths!");}
			break;
		}
		
		
		/**
		 * Tally the read counts each cell barcode receives
		 * Intended to be used before UMIMerging
		 * 		I=... O=...
		 * 
		 * Input arguments
		 * 		I: path to cell barcode corrected file (txt.gz format)
		 * 		O: path to store output file (txt.gz format)
		 * Merge the following cases:
		 * 		Exact matches
		 * 		Match with 1 hamming distance to exactly one other UMI
		 * This means that barcodes that match at 1 hamming distance to more than 1 other sequences are discarded 
		 */
		
		case "ReadcountHistogram":
		{
			// Parse the arguments
			String I = "", O = "";
			for (int i=1; i<args.length; i++)
			{
				String[] j = args[i].split("=");
				
				if (j.length < 2) {throw new java.lang.IllegalArgumentException("Whitespaces are not allowed between argument key specifier and the argument!");}
				
				switch (j[0])
				{
				case "I": I = j[1]; break;
				case "O": O = j[1]; break;
				default: throw new java.lang.IllegalArgumentException("Invalid argument key specifier!");
				}
			}
			
			try (
					BufferedReader brI = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(I)))); // aligned reads
					BufferedWriter bwout = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(O)))) // output file
				)
			{
				
				System.out.printf("%s   ReadcountHistogram   Start%n", LocalDateTime.now().format(time_formatter));
				
				// Create Hash Multiset to get counts of each cell barcode
				Multiset<String> readcounts = HashMultiset.create();
				
				// Read file
				String line;
				while ( (line = brI.readLine()) != null )
				{
					String[] values = line.split(",");
					readcounts.add(values[0]);
				}
				
				// Write to output file, tab-separated
				// First row is header, first column is read count, second column is cell barcode
				bwout.write("Readcounts\tCellbarcode"); bwout.newLine();
				for (String i : Multisets.copyHighestCountFirst(readcounts).elementSet())
				{
					bwout.write(String.format("%d", readcounts.count(i)) + "\t" + i); bwout.newLine();
				}
				
				System.out.printf("%s   ReadcountHistogram   Done%n", LocalDateTime.now().format(time_formatter));
				
			} catch (IOException e) {throw new IllegalArgumentException("Invalid file paths!");}
			
			break;
		}
		
		
		/**
		 * Get digital count from UMI merged file, output a tab-separated txt.gz file
		 * 		I=... O=... CELL_BC_LIST=... SUMMARY=... DUPLICATE_EXPORT=... REMOVE_DUPLICATE=...
		 * Input arguments:
		 * 		I: path to the UMI merged file (txt.gz format)
		 * 		O: path to store the digital count matrix (txt.gz format)
		 * 		CELL_BC_LIST: a list of chosen cell barcodes (from knee plot) (txt format). If NONE (default), export all detected cell barcodes
		 * 		HEADER: whether the CELL_BC_LIST has header to be skipeed (default is false)
		 * 		SUMMARY: path to store the summary file (txt format)
		 * 		DUPLICATE_EXPORT: path to store the list (txt.gz format) of duplicated PLA products across cells (ie, PLA molecules that have the same UMI and PLA IDs).
		 * 				By default, the export file is saved in the current working directory as duplicate_export.txt
		 * 		REMOVE_DUPLICATE: whether to remove duplicated PLA products across cells (default is false)
		 */
		
		case "DigitalCount":
		{
			// Time stamp
			System.out.printf("%s   DigitalCount   Start%n", LocalDateTime.now().format(time_formatter));

			// Parse the arguments
			String I = "", O = "", cell_BC_path = "NONE";
			boolean skip_header = false;
			String SUMMARY = System.getProperty("user.dir") + File.separator + "DigitalCountDropSeq_summary.txt"; // default summary file directory
			String duplicate_export_path = System.getProperty("user.dir") + File.separator + "duplicate_export.txt.gz";
			boolean remove_duplicate = false;
			for (int i=1; i<args.length; i++)
			{
				String[] j = args[i].split("=");
				
				if (j.length < 2) {throw new java.lang.IllegalArgumentException("Whitespaces are not allowed between argument key specifier and the argument!");}
				
				switch (j[0])
				{
				case "I": I = j[1]; break;
				case "O": O = j[1]; break;
				case "CELL_BC_LIST": cell_BC_path = j[1]; break;
				case "HEADER": skip_header = "true".equalsIgnoreCase(j[1]); break;
				case "SUMMARY": SUMMARY = j[1]; break;
				case "DUPLICATE_EXPORT": duplicate_export_path = j[1]; break;
				case "REMOVE_DUPLICATE": remove_duplicate = "true".equalsIgnoreCase(j[1]); break;
				default: throw new java.lang.IllegalArgumentException("Invalid argument key specifier!");
				}
			}
			
			// Read the list of chosen cell barcodes
			// hash set to store chosen cell barcodes, ordered by read counts from drop-seq tools (0 index = most reads)
			Set<String> chosen_BC = new LinkedHashSet<>();
			String BCline;
			
			if (!Objects.equals(cell_BC_path,"NONE"))
			{
				try (
						BufferedReader brBC = new BufferedReader(new FileReader(cell_BC_path)); // cell barcode list
					)
				{
					
					if (skip_header) {brBC.readLine();}
					while ((BCline = brBC.readLine()) != null)
					{
						chosen_BC.add(BCline);
					}
					System.out.printf("%s   DigitalCount   Done processing cell barcode list%n", LocalDateTime.now().format(time_formatter));
					System.out.printf("\tNumber of cell barcodes: %,7d%n", chosen_BC.size());
					
				} catch (IOException e) {throw new IllegalArgumentException("Invalid CELL_BC_LIST argument!");}
			}
			
			// build a list of duplicated PLA products
			System.out.printf("%s   DigitalCount   Building list of duplicated PLA products%n", LocalDateTime.now().format(time_formatter));
			Multiset<String> PLAduplicate = HashMultiset.create(); // count the occurrences of each unique UMI+PLA product across all cell barcodes
			try (
					BufferedReader brI = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(I)))); // aligned reads
					BufferedWriter bwdup = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(duplicate_export_path))))
				)
			{
				String line;
				while ((line = brI.readLine()) != null)
				{
					String[] values = line.split(",");
					if (Objects.equals(cell_BC_path,"NONE") || chosen_BC.contains(values[0]))
					{
						PLAduplicate.add(values[1] + "_" + values[2] + ":" + values[3]); // count occurrences of each UMI+PLA product combination
					}
				}
				
				// Export duplicated PLA products
				bwdup.write("PLA product\tNumber of occurences"); bwdup.newLine();
				for (String i : Multisets.copyHighestCountFirst(PLAduplicate).elementSet())
				{
					if (PLAduplicate.count(i) == 1)
					{
						break;
					}
					bwdup.write(String.format("%s\t%d", i, PLAduplicate.count(i))); bwdup.newLine();
				}
				
			} catch (IOException e) {throw new IllegalArgumentException("Invalid I or DUPLICATE_EXPORT argument!");}
			
			try (
					BufferedReader brI = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(I)))); // aligned reads
					BufferedWriter bwout = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(O)))); // output file
					BufferedWriter bwsum = new BufferedWriter(new FileWriter(SUMMARY)) // summary file
				)
			{
				// Write out the command line arguments
				for (String j : args)
				{
					bwsum.write(j); bwsum.newLine();
					System.out.println(j);
				}
				bwsum.newLine();
				System.out.printf("%n");
				
				
				// Counter for detected duplicated PLA products
				int duplicate_counter = 0;
				
				// Read the PLA products, only keep chosen cell barcodes in CELL_BC_LIST
				String line;
				System.out.printf("%s   DigitalCount   Counting PLA products%n", LocalDateTime.now().format(time_formatter));
				Multiset<String> PLAproduct = HashMultiset.create(); // count the UMIs for each combination of cell barcode-PLA pair
				Set<String> PLA_ID = new TreeSet<>(); // hash set to store unique PLA pairs, ordered by AB1 ID alphabetically
				while ((line = brI.readLine()) != null)
				{
					String[] values = line.split(",");
					if (Objects.equals(cell_BC_path,"NONE") || chosen_BC.contains(values[0]))
					{
						// Add detected cell barcode if CELL_BC_LIST=NONE
						if (Objects.equals(cell_BC_path,"NONE"))
						{
							chosen_BC.add(values[0]);
						}
						
						// Remove duplicate if required
						if (PLAduplicate.count(values[1] + "_" + values[2] + ":" + values[3])> 1)
						{
							duplicate_counter++;
							if (remove_duplicate)
							{continue;}
						}

						PLAproduct.add(values[0] + "_" + values[2] + ":" + values[3]); // count of cell barcode+PLA product combination
						PLA_ID.add(values[2] + ":" + values[3]); // PLA product name

					}
				}
				
				
				System.out.printf("%s   DigitalCount   Writing count matrix%n", LocalDateTime.now().format(time_formatter));
				
				// Write the digital count matrix row by row (ie, by PLA pair)
				// Header: PLA_pair \t cell barcode 1 \t cell barcode 2 ...
				bwout.write("PLA_pair");
				for (String str_i : chosen_BC)
				{
					bwout.write("\t" + str_i);
				}
				bwout.newLine();
				
				// Write each row
				for (String str_i : PLA_ID)
				{
					bwout.write(str_i);
					for (String str_j : chosen_BC)
					{
						bwout.write(String.format("\t%d", PLAproduct.count(str_j + "_" + str_i)));
					}
					bwout.newLine();
				}
				
				// Time stamp
				System.out.printf("%s   DigitalCount   Done%n", LocalDateTime.now().format(time_formatter));
				
				// Write to summary file
				bwsum.write("Number of cell barcodes: " + String.format("%,d", chosen_BC.size())); bwsum.newLine();
				bwsum.write("Number of detected duplicated PLA products: " + String.format("%,d", duplicate_counter)); bwsum.newLine();
				
			} catch (IOException e) {throw new IllegalArgumentException("Invalid file paths!");}
						
			break;
		}
				
		
		/*
		 * Temporary function for debugging
		 * Output reads that have non-matching  connector
		 */
		case "TempDebug":
{
			
			// Parse the arguments
			String R1List = "", ABfile = "", O = "";
			String SUMMARY = System.getProperty("user.dir") + File.separator + "ReadAlignmentSmartSeq_summary.txt"; // default summary file directory
			boolean skip_header = false;
			for (int i=1; i<args.length; i++)
			{
				String[] j = args[i].split("=");
				
				if (j.length < 2) {throw new java.lang.IllegalArgumentException("Whitespaces are not allowed between argument key specifier and the argument!");}
				
				switch (j[0])
				{
				case "R1List": R1List = j[1]; break;
				case "O": O = j[1]; break;
				case "ABfile": ABfile = j[1]; break;
				case "SUMMARY": SUMMARY = j[1]; break;
				case "HEADER": skip_header = "true".equalsIgnoreCase(j[1]); break;
				default: throw new java.lang.IllegalArgumentException("ReadAlignmentSmartSeq: Invalid argument key specifier!");
				}
			}
			
//			System.out.println(R1List);
		
			try (
					BufferedReader br1List = new BufferedReader(new FileReader(R1List)); // List of input read 1 files and the cell barcodes 
					BufferedReader brAB = new BufferedReader(new FileReader(ABfile)); // AB-DNA barcode look up table
					BufferedWriter bwout = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(O)))); // output file
					BufferedWriter bwsum = new BufferedWriter(new FileWriter(SUMMARY)) // summary file
				 )
			{
				
				// Write out the command line arguments
				System.out.println();
				for (String j : args)
				{
					bwsum.write(j); bwsum.newLine();
					System.out.println(j);
				}
				bwsum.newLine();
				System.out.println();
			
				
				// Read the AB look up table into an array
				List<List<String>> ABarray = new ArrayList<List<String>>();
				String ABline;
				if (skip_header) {brAB.readLine();} // skip first line
				while ((ABline=brAB.readLine()) != null)
				{
					String[] values = ABline.split(",");
					ABarray.add(Arrays.asList(values));
				}
	
				
				/**
				 * Process read 1
				 * Read 1 contains PLA products
				 * UMI region: 2nd base to 17th base (16-base long)
				 * Output format: <cell barcode> , <UMI> , <AB1 ID> , <AB2 ID>
				 * Cell barcode is equal to sample barcode
				 */

				// Set up the (total) counters for all read 1 files for the summary file
				
				int read_counter = 0; // read number counter for all fastq files
				int PLA_counter = 0; // PLA product counter
				int bad_connector_counter = 0; // counter for reads with non-matching connector
				int bad_UMI_counter = 0; // counter of reads with bad UMI region
				int non_matching_AB_counter = 0; // counter of the number of reads with non-matching AB barcode
				
				// Set up for alignment
				String R1List_temp; // current read 1 file
				String connector = "TCGTGTCGTGTCGTGTCTAAAG"; // connector sequence
				String[] R1List_temp_array; // array to store each line in R1List
				
				// Start reading
				long my_timer = System.currentTimeMillis();
				System.out.printf("%s   Start alignment%n", LocalDateTime.now().format(time_formatter));
				bwsum.write("Start alignment at " + LocalDateTime.now()); bwsum.newLine(); bwsum.newLine();
				
				while ((R1List_temp=br1List.readLine()) != null)
				{
					R1List_temp_array = R1List_temp.split(",");
					
					
					try (
							BufferedReader br1 = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(R1List_temp_array[0])))); // read1
						)
					{
						// Set up the counters for each read 1 file for the summary file
						int counter = 0; // line number counter in fastq file
//						int PLA_counter_temp = 0; // PLA product counter
//						int bad_connector_counter_temp = 0; // counter for reads with non-matching connector
//						int bad_UMI_counter_temp = 0; // counter of reads with bad UMI region
//						int non_matching_AB_counter_temp = 0; // counter of the number of reads with non-matching AB barcode
						
						String line1; // current line in current read 1 file
						while ((line1=br1.readLine()) != null)
						{
							if ((counter % 4) == 1)
							{			
								
								read_counter++;
								
								// Check if the UMI region has an excessive amount of A
								if (StringUtils.countMatches(line1.substring(1, 17), "A") >= 10)
								{
									if ((read_counter % 1_000_000) == 0)
									{
										System.out.printf("%s   ReadAlignmentSmartSeq   Processed %,15d reads   Elapsed time for last 1,000,000 reads: %ds%n",
												LocalDateTime.now().format(time_formatter), read_counter, (System.currentTimeMillis()-my_timer)/1000);
										my_timer = System.currentTimeMillis();
									}
									
									counter++;
									bad_UMI_counter++;
//									bad_UMI_counter_temp++;
									continue;
								}
								
								int connector_start = 39; // expected starting location of connector is at index 39 (base 40th)
								int connector_start_temp = 39; // temporary starting location of connector, for used in for loop
								
								// Locate the connector using Levenshtein distance, max allowed distance is 2
								// Does not need to take into account N, since the connector region doesn't usually contain N
								int[] connector_shift = new int[] {0, -1, 1}; // only allow up to 1-base shift of the expected starting location
								boolean match_connector = false; // true if matching connector is found
								int connector_distance = 3; // lowest Levenshtein distance found between the true connector and the test connector sequence, will be updated during the for loop
								for (int shift_i : connector_shift)
								{
									int temp_distance = LevenshteinDistance.getDefaultInstance().apply(line1.substring(connector_start+shift_i, connector_start+shift_i+connector.length()), connector);
									if ((temp_distance <= 2) && (temp_distance < connector_distance))
										// match within 2 Levenshtein distance, and only if the temp_distance is lower than previous iterations
									{
										connector_distance = temp_distance;
										connector_start_temp = connector_start + shift_i;
										match_connector = true;
										
										if (temp_distance == 0)
										{ break; }
									}
								}
								connector_start = connector_start_temp;
								
								if (match_connector)
								{
									
									// Check if there is a frameshift, in order to locate AB2 correctly
									// Find the location of TAAAG in the found connector region
									int shift_j = line1.substring(connector_start, connector_start+connector.length()+2).indexOf("TAAAG"); // add 2 just in case there are 2 insertions
									if (shift_j == -1) // skip read if can't find TAAAG in the connector region
									{
										counter++;
										bad_connector_counter++;
										bwout.write(line1);
										bwout.newLine();
//										bad_connector_counter_temp++;
										continue;
									}
									
								}
								else
								{
									bad_connector_counter++;
//									bad_connector_counter_temp++;
									bwout.write(line1);
									bwout.newLine();
								}
								
								
								
								if ((read_counter % 1_000_000) == 0)
								{
									System.out.printf("%s   ReadAlignmentSmartSeq   Processed %,15d reads   Elapsed time for last 1,000,000 reads: %ds%n",
											LocalDateTime.now().format(time_formatter), read_counter, (System.currentTimeMillis()-my_timer)/1000);
									my_timer = System.currentTimeMillis();
								}
							}
							
							
//							if ((((counter-1)/4)+1)>2000000) {System.out.printf("%s   %s   Processed %,d lines%n", LocalDateTime.now().format(time_formatter), counter); break;} // for testing purposes
						
							counter++;
						}
						
//						bwsum.write(R1List_temp_array[0]); bwsum.newLine();
//						bwsum.write("Number of reads with a valid PLA product: " + String.format("%,d", PLA_counter_temp)); bwsum.newLine();
//						bwsum.write("Number of reads discarded because of non-matching connector sequence: " + String.format("%,d",bad_connector_counter_temp)); bwsum.newLine();
//						bwsum.write("Number of reads discarded because of bad UMI: " + String.format("%,d",bad_UMI_counter_temp)); bwsum.newLine();
//						bwsum.write("Number of reads discarded because of non-matching antibody barcode: " + String.format("%,d",non_matching_AB_counter_temp)); bwsum.newLine();
//						bwsum.newLine();
						
					} catch (IOException e) {throw new IllegalArgumentException("Invalid file paths in R1List!");}
					
					
				}
				
				
				System.out.printf("%s   ReadAlignmentSmartSeq   Done: processed %,d reads%n", LocalDateTime.now().format(time_formatter), read_counter);
				System.out.printf("\tNumber of reads with a valid PLA product: %,15d%n", PLA_counter);
				System.out.println();
				
				// Write to summary file
				bwsum.write("TempDebug: Finished at " + LocalDateTime.now().withNano(0) + ", processed " + String.format("%,d",read_counter) + " reads"); bwsum.newLine();
				bwsum.write("Total number of reads with a valid PLA product: " + String.format("%,d", PLA_counter)); bwsum.newLine();
				bwsum.write("Total number of reads discarded because of non-matching connector sequence: " + String.format("%,d",bad_connector_counter)); bwsum.newLine();
				bwsum.write("Total number of reads discarded because of bad UMI: " + String.format("%,d",bad_UMI_counter)); bwsum.newLine();
				bwsum.write("Total number of reads discarded because of non-matching antibody barcode: " + String.format("%,d",non_matching_AB_counter)); bwsum.newLine();
				
			} catch (IOException e) {throw new IllegalArgumentException("Invalid file paths!");}
			break;
		}
		
		default: throw new java.lang.IllegalArgumentException("Invalid function call!");
		}
	}

}