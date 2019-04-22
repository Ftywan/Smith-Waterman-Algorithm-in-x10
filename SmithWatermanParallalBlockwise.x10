package smithwatermanalgo;
import x10.lang.String;
import x10.io.Console;
import x10.io.File;
import x10.io.EOFException;
import x10.util.StringBuilder;
import x10.io.FileReader;
import x10.lang.Exception;
import x10.util.HashMap;
import x10.array.Array_2;
import x10.lang.Math;
import x10.io.Console;
import x10.lang.Char;
import x10.lang.*;
import x10.io.*;
import x10.xrx.Runtime;



public class SmithWatermanParallalBlockwise {

    private val NUM_COLS_IN_BLOCK:Int = 16n;
    private val NUM_ROWS_IN_BLOCK:Int;

    private val NUM_BLOCKS_X:Int = Runtime.NTHREADS;
    private val NUM_BLOCKS_Y:Int;

    private var finishStatus:Array_2[Int];

    //two sequences of AA
    private val seq1:String;
    private val seq2:String;
    private val length1:Int;
    private val length2:Int;

    private var maxi:Int;
    private var maxj:Int;

    private var outstr1arr:Rail[Char];
    private var outstr2arr:Rail[Char];
    private var lengthOut:Int;

    //the score matrix
    private var score:Array_2[Int];

    private var scoreLeft:Array_2[Int];
    private var scoreUp:Array_2[Int];

    //similarity fuction constant
    private val SCORE_THRESHOLD:Int = 19n;
    private val GAP_OPENING_PANALTY:Int;
    private val GAP_EXTENSION_PANALTY:Int;
    private val INDEL_SCORE:Int = -9n;

    //constant in direction matrix
    private val DR_LEFT:Int = 1n;// 0001
    private val DR_UP:Int = 2n;// 0010
    private val DR_DIAG:Int = 4n;// 0100
    private val DR_ZERO:Int = 8n;// 1000



    //direction matrix
    private var prevCells:Array_2[Int];

    private val blosum62:Array_2[Int];

    private val seqToNum:HashMap[Char, Int];
    private val blosumReader:BlosumReader;
    private val fastaReader:FastaReader;
    

    public def this(fastaName1:String, fastaName2:String, blosumFileName:String, openPanalty:Int, extensionPanalty:Int) {
        
        fastaReader = new FastaReader();
        this.blosumReader = new BlosumReader(blosumFileName);
        this.blosum62 = blosumReader.getBlosum62();
        this.seqToNum = blosumReader.getSeqToNum();


        this.seq1 = fastaReader.readFastaFile(fastaName1);
        this.seq2 = fastaReader.readFastaFile(fastaName2);
        this.length1 = seq1.length();
        this.length2 = seq2.length();

        Console.OUT.println(length1);

        this.GAP_OPENING_PANALTY = openPanalty;
        this.GAP_EXTENSION_PANALTY = extensionPanalty;

        score = new Array_2[Int](length1 + 1, length2 + 1);
        scoreLeft = new Array_2[Int](length1 + 1, length2 + 1);
        scoreUp = new Array_2[Int](length1 + 1, length2 + 1);
        prevCells = new Array_2[Int](length1 + 1, length2 + 1);

        this.NUM_ROWS_IN_BLOCK = Math.ceil((this.length1 as Double) / this.NUM_BLOCKS_X) as Int;
        this.NUM_BLOCKS_Y = Math.ceil((this.length2 as Double) / this.NUM_COLS_IN_BLOCK) as Int;

        this.finishStatus = new Array_2[Int](NUM_BLOCKS_X, NUM_BLOCKS_Y);

        for (i in 0..(NUM_BLOCKS_X - 1n)) {
            for (j in 0..(NUM_BLOCKS_Y -1n )) {
                finishStatus(i, j) = 0n;
            }
        }
    }

    public def similarity(i:Int, j:Int):Int {
        //return blosum62(seqToNum.get(seq1.charAt(i-1n)), seqToNum.get(seq2.charAt(j-1n)));
        var char1:Char = seq1.charAt(i-1n);
        var char2:Char = seq2.charAt(j-1n);
        if (seqToNum.containsKey(char1)) {
            if (seqToNum.containsKey(char2)) {
                return blosum62(seqToNum.get(char1), seqToNum.get(char2));
            } else {
                return blosum62(seqToNum.get(char1), 23n);
            }
        } else {
            if (seqToNum.containsKey(char2)) {
                return blosum62(23n, seqToNum.get(char2));
            } else {
                return blosum62(23n, 23n);
            }
        }
    }

    public def buildMatrix() {

        for (i in 0..(NUM_BLOCKS_X - 1n)) {
            for (j in 0..(NUM_BLOCKS_Y -1n )) {
                var point:Rail[Int] = getBlockPosition(i as Int, j as Int);
                Console.OUT.println("i:" + i + " j:" + j + " " + point(0n) + " " + point(1n) + " " + point(2n) + " " + point(3n));
            }
        }

        var max:Int = -99999999n;
        var maxi:Int = -1n;
        var maxj:Int = -1n;

        //base case
        score(0, 0) = 0n;
        scoreLeft(0, 0) = -9999990n;
        scoreLeft(0, 0) = -9999990n;
        prevCells(0,0) = DR_ZERO;

        //the first row
        for(i in 1..length1) {
            score(i, 0) = 0n;
            scoreLeft(i, 0) = -9999990n;
            scoreUp(i, 0) = -9999990n;
            prevCells(i, 0) = DR_ZERO;
        }

        //the first column
        for(j in 1..length2) {
            score(0, j) = 0n;
            scoreLeft(0, j) = -9999990n;
            scoreUp(0, j) = -9999990n;
            prevCells(0, j) = DR_ZERO;
        }

        finish for (i in 0n..(this.NUM_BLOCKS_X - 1n)) async {
                var maxResult:Rail[Int] = workerThread(i);
                Console.OUT.println(i + " returned.");
                if (maxResult(2n) > max) {
                    max = maxResult(2n);
                    maxi = maxResult(0n);
                    maxj = maxResult(1n);
                }
        }

        //var point:Rail[Int] = diagnalCover(1n, 1n, length1, length2);

        this.maxi = maxi;
        this.maxj = maxj;
    }

    public def workerThread(var id:Int):Rail[Int] {
        Console.OUT.println(id + " started.");
        var max:Int = -99999999n;
        var maxi:Int = -1n;
        var maxj:Int = -1n;
        if (id == 0n) {
            for (i in 0n..(this.NUM_BLOCKS_Y - 1n)) {
                var points:Rail[Int] = getBlockPosition(id, i);
                var maxResult:Rail[Int] = diagnalCover(points(0n), points(1n), points(2n), points(3n));
                if (maxResult(2n) > max) {
                    max = maxResult(2n);
                    maxi = maxResult(0n);
                    maxj = maxResult(1n);
                }
                this.finishStatus(id, i) = 1n;
                Console.OUT.println(id + " " + i + " finished");
            }
        } else {
            Console.OUT.println(id + " ready");
            for (i in 0n..(this.NUM_BLOCKS_Y - 1n)) {
                Console.OUT.println(id + " " + i + " waiting");
                while (this.finishStatus(id -1n, i) != 1n) {}
                var points:Rail[Int] = getBlockPosition(id, i);
                var maxResult:Rail[Int] = diagnalCover(points(0n), points(1n), points(2n), points(3n));
                if (maxResult(2n) > max) {
                    max = maxResult(2n);
                    maxi = maxResult(0n);
                    maxj = maxResult(1n);
                }
                this.finishStatus(id, i) = 1n;
                Console.OUT.println(id + " " + i + " finished");
            }
        }
        var point:Rail[Int] = [maxi, maxj, max];
        return point;
    }

    def getBlockPosition(var i:Int, var j:Int):Rail[Int] {

        var position:Rail[Int];

        //divide the score matrix to blocks
        //get the num of blocks
        //size of the martix is (length+1n)(length2+1n)
        /*
        var numOfBlocksInHight = length1 / NUM_ROWS_IN_BLOCK + 1n;
        var numOfBlocksInWidth = length2 / NUM_COLS_IN_BLOCK + 1n;
        */

        //get the final position
        var leftI:Int;
        var leftJ:Int;
        var rightI:Int;
        var rightJ:Int;

        //topleft position: never excel, need not consider the special case
        leftI = i * NUM_ROWS_IN_BLOCK + 1n;
        leftJ = j * NUM_COLS_IN_BLOCK + 1n;

        //special case: edge of the matrix
        if(i == NUM_BLOCKS_X-1n) {
            rightI = length1;
        }else {
            rightI = (i + 1n) * NUM_ROWS_IN_BLOCK + 1n;
        }

        if(j == NUM_BLOCKS_Y-1n) {
            rightJ = length2;
        }else {
            rightJ = (j + 1n) * NUM_COLS_IN_BLOCK + 1n;
        }

        position = [leftI, leftJ, rightI, rightJ];

        return position;
    }

    public def diagnalCover(var a1:Int, var b1:Int, var a2:Int, var b2:Int):Rail[Int] {
        var max:Int = -99999999n;
        var maxi:Int = -1n;
        var maxj:Int = -1n;
        var temp:Int = -1n;
        //start from row 0
        for(j in b1 .. b2) {
            var i:Int = a1;
            var k:Int = j;
            while(i <= a2 && k >= b1) {
                temp = calculateScore(i, k);
                if (temp > max) {
                    max = temp;
                    maxi = i;
                    maxj = k;
                }
                i++;
                k--;
            }
        }

        //continue from final col
        for(i in (a1+1n) .. a2) {
            var j:Int = b2;
            var k:Int = i;
            while(k <=a2 && j >= b1) {
                temp = calculateScore(k, j);
                if (temp > max) {
                    max = temp;
                    maxi = k;
                    maxj = j;
                }
                k++;
                j--;
            }
        }

        var point:Rail[Int] = [maxi, maxj, max];

        return point;
    }


    public def calculateScore(var i:Int, var j:Int):Int {
        //Console.OUT.println(i + " " + j);
        var diagScore:Int = score(i-1, j-1) + similarity(i, j);
        //Console.OUT.println("dig:" + score(i-1, j-1) + " simi:" + similarity(i, j));

        var newOpenGapLeftScore:Int = score(i, j-1) - GAP_OPENING_PANALTY;
        var newExtentionGapLeftScore:Int = scoreLeft(i, j-1) - GAP_EXTENSION_PANALTY;
        scoreLeft(i, j) = Math.max(newOpenGapLeftScore, newExtentionGapLeftScore);

        var newOpenGapUpScore:Int = score(i-1, j) - GAP_OPENING_PANALTY;
        var newExtentionGapUpScore:Int = scoreUp(i-1, j) - GAP_EXTENSION_PANALTY;
        scoreUp(i, j) = Math.max(newOpenGapUpScore,newExtentionGapUpScore);


        var upScore:Int = scoreUp(i, j);
        var leftScore:Int = scoreLeft(i, j);

        //score(i, j) = Math.max(diagScore, Math.max(upScore, Math.max(leftScore, 0n)));
        score(i, j) = Math.max(diagScore, Math.max(upScore, Math.max(leftScore, 0n)));
        prevCells(i, j) = 0n;

        if (diagScore == score(i, j)) {
            prevCells(i, j) |= DR_DIAG;
        } else if (leftScore == score(i, j)) {
            prevCells(i, j) |= DR_LEFT;
        } else if (upScore == score(i, j)) {
            prevCells(i, j) |= DR_UP;
        } else if (0n == score(i, j)) {
            prevCells(i, j) |= DR_ZERO;
        }
        //Console.OUT.println("final:" + score(i, j));
        return score(i, j);
    }

    public def getMaxScore():Int {
        var maxScore:Int = this.score(maxi, maxj);
        /*
        for(i in 1 .. length1) {
            for(j in 1 .. length2) {
                if(score(i, j) > maxScore) {
                    maxScore = score(i, j);
                }
            }
        }
        */
        return maxScore;
    }


    // TODO: printAlignments()

    //returns the end point of tracing back (the top left cell) and the number of matches
    public def traceback(var i:Int, var j:Int):Rail[Int] {
        var num:Int = 0n;
        var match:Int = 0n;
        var identity:Int = 0n;
        var gap:Int = 0n;
        var traceSTR1:Rail[Char] = new Rail[Char](Math.max(length1, length2)*2);
        var traceSTR2:Rail[Char] = new Rail[Char](Math.max(length1, length2)*2);
        var str1len:Int = 0n;
        var str2len:Int = 0n;
        this.outstr1arr = new Rail[Char](Math.max(length1, length2)*2);
        this.outstr2arr = new Rail[Char](Math.max(length1, length2)*2);
        //find the direction to traceback
        while (true)
        {
            if ((prevCells(i, j) & DR_LEFT) > 0n) {
                num ++;
                gap ++;
                traceSTR1(str1len) = '-';
                str1len ++;
                traceSTR2(str2len) = seq2.charAt(j-1n);
                str2len ++;
                j--;
                //if (score(i-1n, j)>0n) i--;
                //else    break;              
            }
            if ((prevCells(i, j) & DR_UP) > 0n) {
                num ++;
                gap ++;
                traceSTR1(str1len) = seq1.charAt(i-1n);
                str1len ++;
                traceSTR2(str2len) = '-';
                str2len ++;
                i--;
//          return traceback(i, j-1);
                //if (score(i, j-1n)>0n) j--;
                //else    break;              
            }
            if ((prevCells(i, j) & DR_DIAG) > 0n) {
                if (seq1.charAt(i-1n) == seq2.charAt(j-1n)) {
                    identity ++;
                }
                num ++;
                match ++;
                traceSTR1(str1len) = seq1.charAt(i-1n);
                traceSTR2(str2len) = seq2.charAt(j-1n);
                str1len ++;
                str2len ++;
                j--;
                i--;
//          return traceback(i-1, j-1);
                //if (score(i-1n, j-1n)>0n) {i--;j--;}
                //else     break;             
            } else {
                break;
            }
        }

        for (q in 0..(str1len-1)) {
            this.outstr1arr(str1len - q - 1) = traceSTR1(q);
            this.outstr2arr(str1len - q - 1) = traceSTR2(q);
        }

        this.lengthOut = str1len;
        var point:Rail[Int] = [i, j, num, match, gap, identity];
        return point;
    }

    public def printStr1() {
        for (i in 0..(lengthOut-1)) {
            Console.OUT.print(this.outstr1arr(i));
        }
    }

    public def printStr2() {
        for (i in 0..(lengthOut-1)) {
            Console.OUT.print(this.outstr2arr(i));
        }
    }

    public def getMatchNumber():Int {

        var matches:Int = 0n;

        for(i in 1n..length1) {
            for(j in 1n..length2) {
                if(score(i, j) > SCORE_THRESHOLD && score(i, j) > score(i-1n, j)
                    && score(i, j) > score(i ,j-1n) && score(i, j) > score(i-1n, j-1n))
                {
                    if(i == length1 || j == length2 || score(i, j) > score(i+1n, j+1n))
                    {
                        var endPoint:Rail[Int] = traceback(i, j);

                        matches += endPoint(2n);
                    }
                }
            }
        }
        return matches;
    }

    public def getMatchGap():Rail[Int] {
        /*
        var max:Int = 0n;
        var maxJ:Int = 0n;
        for(j in 1n..length2) {
            if (score(length1, j) > max) {
                max = score(length1, j);
                maxJ = j;
            } 
        }
        */
        //var endPoint:Rail[Int] = traceback(length1, maxJ);
        var endPoint:Rail[Int] = traceback(this.maxi, this.maxj);
        var result:Rail[Int] = [endPoint(3n), endPoint(4n), endPoint(5n), endPoint(2n)];
        return result;
    }

    public def getTototalNumber():Int {
        return 0n;
    }

    public def getGapNumber():Int {
        return 0n;
    }

    public def printAlignments() {

    }

    public def gettotalScore():Int {
        return 0n;
    }

    public def printMatrix() {
        for (i in 1..length1) {
            for (j in 1..length2) {
                Console.OUT.print("|(" + i + " " + j +") " +score(i, j) + "| ");
            }
            Console.OUT.print("\n");
        }
        Console.OUT.print("\n");
        for (i in 1..length1) {
            for (j in 1..length2) {
                Console.OUT.print("|(" + i + " " + j +") " +scoreLeft(i, j) + "| ");
            }
            Console.OUT.print("\n");
        }
        Console.OUT.print("\n");
        for (i in 1..length1) {
            for (j in 1..length2) {
                Console.OUT.print("|(" + i + " " + j +") " +scoreUp(i, j) + "| ");
            }
            Console.OUT.print("\n");
        }
    }


    public static def main(argv: Rail[String]) {
        Console.OUT.println("Input the FASTA_FILE_1 FASTA_FILE_2 MATCH_FILE GAP_OPENING_PANALTY GAP_EXTENSION_PANALTY");
        //Console.OUT.println(Runtime.NTHREADS);
        //Console.OUT.println(Runtime.MAX_THREADS);
        val s = x10.io.Console.IN.readLine();

        val param = s.split(" ");

        val fasta1:String = param(0n);
        val fasta2:String = param(1n);
        val match:String = param(2n);
        val openPanalty:Int = Int.parse(param(3n));
        val extPanalty:Int = Int.parse(param(4n));


        val sw:SmithWatermanParallalBlockwise = new SmithWatermanParallalBlockwise(fasta1, fasta2, match, openPanalty, extPanalty);
        sw.buildMatrix();
        //sw.printMatrix();

        //Console.OUT.println("IO debug");
        //Console.OUT.println(sw.seq1);
        //Console.OUT.println(sw.seq2);
        //Console.OUT.println(sw.blosumFileName);
        //Console.OUT.println(sw.GAP_OPENING_PANALTY);
        //Console.OUT.println(sw.GAP_EXTENSION_PANALTY);
        var result:Rail[Int] = sw.getMatchGap();

        Console.OUT.println("Identity: " + result(2) + "/" + result(3) + " (" + ((result(2) as Double)/ result(3)) + ")");

        Console.OUT.println("Gaps: " + result(1) + "/" + result(3) + " (" + ((result(1) as Double) / result(3)) + ")");

        Console.OUT.println("Score: " + sw.getMaxScore());

        Console.OUT.print("1: ");
        sw.printStr1();
        Console.OUT.println("\n");

        Console.OUT.print("2: ");
        sw.printStr2();
        Console.OUT.println("\n");

    }

}







class FastaReader {
    public static def readFastaFile(fastaFileName:String):String {
        val fastaFile:File = new File(fastaFileName);
        val fastaReader:FileReader = fastaFile.openRead();
        val header:String = fastaReader.readLine();
        val builder:StringBuilder = new StringBuilder();
        var line:String = null;
        var iterator:ReaderIterator[String] = fastaReader.lines();
        while (iterator.hasNext()) {
            line = iterator.next();
            //Console.OUT.println(line);
            line = line.trim();
            builder.add(line);
        }
        return builder.toString();
    } 
}

class BlosumReader {
    public var BLOSUM62: Array_2[Int];
    private var SeqToNum: HashMap[Char, Int];
    private var NumToSeq: HashMap[Int, Char];
    private var NUMOFSEQ: Int = 23n;

    public def this(blosumFileName: String) {
        this.BLOSUM62 = new Array_2[Int](NUMOFSEQ + 1n, NUMOFSEQ + 1n);
        this.readBlosumFile(blosumFileName);
    }


    private def readBlosumFile(blosumFileName: String) {
        val fastaFile: File = new File(blosumFileName);
        val fastaReader: FileReader = fastaFile.openRead();
        val header: String = fastaReader.readLine().trim();
        this.initSeqToNum(header);
        var line: String = null;
        var chars: Rail[String];
        for (i in 0n .. (NUMOFSEQ)) {
            try {
                line = fastaReader.readLine().trim();
           
                
            } catch (e: EOFException) {
                break;
            }
            chars = line.substring(1n).trim().split(" ");
            
            var pointer: Int = 0n;
            for (j in 0n .. (chars.size - 1n)) {
                if (chars(j).equals(" ") || chars(j).length() <= 0) {
                    continue;
                } else {
                    this.BLOSUM62(i, pointer) = Int.parseInt(chars(j));
                   
                    pointer++;
                }
            }
        }
        
    }

    private def initSeqToNum(header:String) {
        this.SeqToNum = new HashMap[Char, Int]();
        this.NumToSeq = new HashMap[Int, Char]();
        var chars: Rail[String] = header.split(" ");
        var value: Int = 0n;
        for (i in (0n .. (chars.size - 1n))) {
            if (chars(i).equals("")) {
                continue;
            } else {
                SeqToNum.put(chars(i).charAt(0n), value);
                NumToSeq.put(value, chars(i).charAt(0n));
            }
            value++;
        }
        this.NUMOFSEQ = value;
        
    }

    public def getSeqToNum(): HashMap[Char, Int] {
        return this.SeqToNum;
    }

    public def getNumToSeq(): HashMap[Int, Char] {
        return this.NumToSeq;
    }

    public def getBlosum62(): Array_2[Int] {
        return this.BLOSUM62;
    }

}







