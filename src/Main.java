import java.io.*;
import java.util.Scanner;
import java.util.ArrayList;

public class Main {
    // 1.1
    //je m'en criss un peu ici so j'ai tout criss dans un tableau, dans le fond
    //les index 0,2,4,6... sont les >dsad|dasd|sdads| et les index 1,3,5,... sont
    //les nucleotides, good luck
    //trop de trouble de les filtrer

    //2e argument la dabatase quon creer dans la main
    public static ArrayList<DBentry> readFasta(String file) throws FileNotFoundException {
        ArrayList<DBentry> database =  new ArrayList<DBentry>();
        try (Scanner sc = new Scanner(new File(file))) {

            while (sc.hasNextLine()) {
                String infos = sc.nextLine();
                String sequence = sc.nextLine();

                //create entry that contains these 2 strings
                DBentry currEntry = new DBentry(infos,sequence);

                //add entry to database
                database.add(currEntry);
            }
        }
        return database;
    }

    //1.2 and 1.3
    //la position de chaque KMER est l'index dans le tableau
    public static ArrayList<String> kmer (String word, int k){
        ArrayList<String> mer = new ArrayList<String>();

        for (int i  = 0; i <= word.length() - k ; i++){
            mer.add(word.substring(i, i+k));
        }

        return mer;
    }

    //retourne une list avec les hsps
    public static ArrayList<ArrayList<String>> hsp(String input, ArrayList<DBentry> data, String seed, ArrayList<ArrayList<Integer>> position, ArrayList<ArrayList<Integer>> positioninput){
        int k = seed.length();
        ArrayList<ArrayList<String>> hsps = new ArrayList<ArrayList<String>>();
        ArrayList<String> kinput = kmer(input, k);

        for (int y = 0; y < data.size(); y++) {
            ArrayList<String> kdata = kmer(data.get(y).getSequence(), k);
            hsps.add(new ArrayList<String>());
            position.add(new ArrayList<Integer>());
            positioninput.add(new ArrayList<Integer>());
            for (int i = 0; i < kinput.size(); i++) {
                for (int j = 0; j < kdata.size(); j++) {
                    Boolean bruh = true;
                    for (int z = 0; z < k; z++) {
                        if (kinput.get(i).charAt(z) != kdata.get(j).charAt(z) && seed.charAt(z) == '1') {
                            bruh = false;
                            break;
                        } else if (kinput.get(i).charAt(z) == kdata.get(j).charAt(z) && seed.charAt(z) == '0') {
                            bruh = false;
                            break;
                        }
                    }
                    if (bruh) {
                        hsps.get(y).add(kinput.get(i));
                        position.get(y).add(j); // y ieme sequence dans la database , j ieme kmer = index de cette sequence
                        positioninput.get(y).add(i);// y ieme sequence dans la database, i eme kmer = index de notre sequence inconnue
                    }
                }
            }
        }
        return hsps;
    }

    //1.3
    //1.4
    public static ArrayList<ArrayList<String>> extendHsps (ArrayList<ArrayList<String>> hsps,ArrayList<DBentry> database, String sequence, ArrayList<ArrayList<Integer>> dbPositions, ArrayList<ArrayList<Integer>> sequencePositions, int seuil) {

        for (int i = 0; i < hsps.size(); i++) { // kmers in this database line

            boolean endLeft = false;
            boolean endRight = false;

            for (int j = 0; j < hsps.get(i).size(); j++) { // current kmer

                //database position
                int kmerDBstartPosition = dbPositions.get(i).get(j);
                int kmerDBendPosition = kmerDBstartPosition + hsps.get(i).get(j).length() - 1;


                //sequence position
                int kmerSequenceStartPosition = sequencePositions.get(i).get(j);
                int kmerSequenceEndPosition = kmerSequenceStartPosition + hsps.get(i).get(j).length() - 1;

                int maxScore = 0;
                int currentScore = 0;


                boolean leftMatch = false;
                boolean rightMatch = false;

                while (currentScore > maxScore - seuil && (!endLeft && !endRight)) {



                    leftMatch = false;
                    rightMatch = false;

                    //outOfBoundsVerification
                    if (kmerDBstartPosition - 1 == -1 || kmerSequenceStartPosition - 1 == -1) {
                        endLeft = true;
                    }
                    if (kmerDBendPosition + 1 == database.get(i).getSequence().length() || kmerSequenceEndPosition + 1 == sequence.length()) {
                        endRight = true;
                    }


                    //left match
                    if (!endLeft && database.get(i).getSequence().charAt(kmerDBstartPosition - 1) == sequence.charAt(kmerSequenceStartPosition - 1)) {
                        leftMatch = true;
                    }
                    //right match
                    if (!endRight && database.get(i).getSequence().charAt(kmerDBendPosition + 1) == sequence.charAt(kmerSequenceEndPosition + 1)) {
                        rightMatch = true;
                    }


                    //left match only
                    if (!endLeft && leftMatch && !rightMatch) {

                        char leftCar = sequence.charAt(kmerSequenceStartPosition -1);

                        //extend hsp
                        hsps.get(i).set(j, leftCar + hsps.get(i).get(j));
                        //extend start positions
                        dbPositions.get(i).set(j, kmerDBstartPosition - 1);
                        sequencePositions.get(i).set(j, kmerSequenceStartPosition - 1);

                        //score acc
                        currentScore += 5;
                        if (currentScore > maxScore) {
                            maxScore = currentScore;
                        }

                        //next car
                        kmerDBstartPosition--;
                        kmerSequenceStartPosition--;
                    }


                    //right match only
                    if (!endRight && rightMatch && !leftMatch) {

                        char rightCar = sequence.charAt(kmerSequenceEndPosition + 1);
                        //extend hsp
                        hsps.get(i).set(j, hsps.get(i).get(j) + rightCar);

                        //score acc
                        currentScore += 5;
                        if (currentScore > maxScore) {
                            maxScore = currentScore;
                        }

                        //next car
                        kmerDBendPosition++;
                        kmerSequenceEndPosition++;
                    }



                    //extend both ways
                    if (!endLeft && !endRight && leftMatch && rightMatch) {

                        char leftCar = sequence.charAt(kmerSequenceStartPosition -1);
                        char rightCar = sequence.charAt(kmerSequenceEndPosition +1);

                        //extend hsp
                        hsps.get(i).set(j, leftCar + hsps.get(i).get(j) + rightCar);

                        //score acc
                        currentScore += 10;

                        if (currentScore > maxScore) {
                            maxScore = currentScore;
                        }

                        //extend start positions
                        dbPositions.get(i).set(j, kmerDBstartPosition - 1);
                        sequencePositions.get(i).set(j, kmerSequenceStartPosition - 1);

                        //next car
                        kmerDBstartPosition--;
                        kmerSequenceStartPosition--;

                        kmerDBendPosition++;
                        kmerSequenceEndPosition++;

                    }

                    //no match -> go left
                    if (!endLeft && !rightMatch && !leftMatch) {


                        char leftCar = sequence.charAt(kmerSequenceStartPosition -1);

                        //extend hsp
                        hsps.get(i).set(j, leftCar + hsps.get(i).get(j));
                        //extend start positions
                        dbPositions.get(i).set(j, kmerDBstartPosition - 1);
                        sequencePositions.get(i).set(j, kmerSequenceStartPosition - 1);


                        //score acc
                        currentScore -= 4;

                        //next car
                        kmerDBstartPosition--;
                        kmerSequenceStartPosition--;
                    }



                }


            }
        }
        return hsps;
    }




    //extend de stalin

    public static <T> void removeDuplicates(ArrayList<ArrayList<T>> list) {

        // Create a new ArrayList
        ArrayList<ArrayList<T>> newList = new ArrayList<ArrayList<T>>();

        // Traverse through the first list
        for (int i = 0;  i < list.size(); i++) {
            if (!list.get(i).isEmpty()) {
                ArrayList<T> newist = new ArrayList<T>();
                for (int j = 0; j < list.get(i).size(); j++) {
                    if (!newist.contains(list.get(i).get(j))) {
                        newist.add(list.get(i).get(j));
                    }
                }
                list.set(i, newist);
            }
        }
    }

    //1.4
    public static ArrayList<ArrayList<String>> glouton (ArrayList<DBentry> database, String input, ArrayList<ArrayList<Integer>> positionInput, ArrayList<ArrayList<Integer>> positionDatabase, ArrayList<ArrayList<String>> hsp, Integer seuil, String seed){
        ArrayList<ArrayList<String>> hsps = hsp;

        for (int i =0; i < hsps.size(); i++){
            if (!hsp.get(i).isEmpty()) {
                for (int j = 0; j < hsp.get(i).size(); j++) {

                    int seqInputFirst = positionInput.get(i).get(j)-1;
                    int seqInputLast = positionInput.get(i).get(j) + seed.length();
                    int seqDataFirst = positionDatabase.get(i).get(j)-1;
                    int seqDataLast = positionDatabase.get(i).get(j) + seed.length();

                    int maxScore = 0;
                    int score = 0;

                    boolean extendLeft = true;
                    boolean extendRight = true;

                    while(extendLeft || extendRight){
//                        System.out.println("-------------------------");
//                        System.out.println(maxScore);
//                        System.out.println(score);
//                        System.out.println(seuil);
                        int scoreLeft = 0;
                        int scoreRight = 0;

                        if (seqInputFirst == -1 || seqDataFirst == -1){
                            extendLeft = false;
                        } else {
                            if (input.charAt(seqInputFirst) == database.get(i).getSequence().charAt(seqDataFirst)){
                                scoreLeft += 5;
                            } else {
                                scoreLeft -= 4;
                            }
                        }

                        if (seqDataLast == database.get(i).getSequence().length()  || seqInputLast == input.length()) {
                            extendRight = false;
                        } else {
                            if (input.charAt(seqInputLast) == database.get(i).getSequence().charAt(seqDataLast)){
                                scoreRight += 5;
                            } else {
                                scoreRight -= 4;
                            }
                        }

                        if (!extendLeft && extendRight) {
                            score += scoreRight;
                            if (maxScore - score < seuil){
                                hsps.get(i).set(j, hsps.get(i).get(j)+input.charAt(seqInputLast));
                                seqInputLast++;
                                seqDataLast++;
                            }
                        } else if (extendLeft && !extendRight) {
                            score += scoreLeft;
                            if (maxScore - score < seuil) {
                                hsps.get(i).set(j, input.charAt(seqInputFirst) + hsps.get(i).get(j));
                                positionDatabase.get(i).set(j, positionDatabase.get(i).get(j) - 1);
                                positionInput.get(i).set(j, positionInput.get(i).get(j) - 1);
                                seqInputFirst--;
                                seqDataFirst--;
                            }
                        } else if (extendLeft && extendRight) {
                            if (scoreLeft >= scoreRight){
                                score += scoreLeft;
                                if (maxScore - score < seuil) {
                                    hsps.get(i).set(j, input.charAt(seqInputFirst) + hsps.get(i).get(j));
                                    positionDatabase.get(i).set(j, positionDatabase.get(i).get(j) - 1);
                                    positionInput.get(i).set(j, positionInput.get(i).get(j) - 1);
                                    seqInputFirst--;
                                    seqDataFirst--;
                                }
                            } else {
                                score += scoreRight;
                                if (maxScore - score < seuil) {
                                    hsps.get(i).set(j, hsps.get(i).get(j) + input.charAt(seqInputLast));
                                    seqInputLast++;
                                    seqDataLast++;
                                }
                            }
                        }

                        if (score > maxScore){
                            maxScore = score;
                        }

                        if (maxScore - score >= seuil){
                            extendLeft = false;
                            extendRight = false;
                        }
                    }
                }
            }
        }

        removeDuplicates(hsps);
        removeDuplicates(positionInput);
        removeDuplicates(positionDatabase);

        return hsps;
    }




    //1.5
    public static ArrayList<ArrayList<String>> fusionHsp (ArrayList<ArrayList<String>> extendedHsps,ArrayList<DBentry> database, ArrayList<ArrayList<Integer>> dbPositions, ArrayList<ArrayList<Integer>> sequencePositions) {

        //condition de fusion : end index dun hsp > start position dun autre && < end


            for (int i = 0; i < extendedHsps.size(); i++) {
                for (int j = 0; j < extendedHsps.get(i).size(); j++) {

                    //sortir si ya aucun hsp a cette ligne ou seulement 1 (donc pas de chevauchement->skip line)
                    if(extendedHsps.get(i).size()==0 || extendedHsps.get(i).size()==1){
                        break;
                    }

                    //end index of current hsp we want to compare
                    int hspSize = extendedHsps.get(i).get(j).length();
                    int hspEndIndex = sequencePositions.get(i).get(j) + hspSize - 1;

                    //chek all other hsps
                    for (int l = 0; l < extendedHsps.get(i).size(); l++) {


//                        System.out.println("i = "+i);
//                        System.out.println("j = " +j);
//                        System.out.println("l = "+l);

                            //start and end index of the other hsp
                            int otherHspStart = sequencePositions.get(i).get(l);
                            int otherHspSize = extendedHsps.get(i).get(l).length();
                            int otherHspEnd = sequencePositions.get(i).get(l) + otherHspSize - 1;

                            //condition de fusion
                            if (hspEndIndex > otherHspStart && hspEndIndex < otherHspEnd) {

                                //a partir de end index du current hsp, concat lautre hsp a part sa partie chevauchee
                                int longueurChevauchement = hspEndIndex - otherHspStart;

                                String concatPartOfOtherHsp = (extendedHsps.get(i).get(l).substring(longueurChevauchement));

                                String fusionFinale = extendedHsps.get(i).get(j) + concatPartOfOtherHsp;


                                //add fusion
                                extendedHsps.get(i).set(j, fusionFinale);

                            }

                    }

                }
            }
        return extendedHsps;
    }






    public static int calculateBruteScore(ArrayList<ArrayList<String>> hsps,ArrayList<DBentry> database,ArrayList<ArrayList<Integer>>databasePositions, int i, int j){

        String currentHsp = hsps.get(i).get(j);
        int dbSequenceNumber = i;
        int sequenceIndex = databasePositions.get(i).get(j);

        int score = 0;
        for (int k = 0; k<currentHsp.length(); k++){
            if (currentHsp.charAt(k) == database.get(dbSequenceNumber).getSequence().charAt(sequenceIndex)){
                score +=5;
            }else{
                score-=4;
            }
            j++;
        }
        return score;

    }

    public static int databaseTotalSize(ArrayList<DBentry> database){
        int len = 0;
        for (int i= 0; i<database.size(); i++){
            len += database.get(i).sequence.length();
        }
        return len;
    }




    //1.6
    public static ArrayList<HspStat> getHspStat(ArrayList<ArrayList<String>> hsps,ArrayList<DBentry> database,String sequence, ArrayList<ArrayList<Integer>>databasePositions,ArrayList<ArrayList<Integer>>sequencePositions ){

        double lamda = 0.192;
        double k = 0.176;

        ArrayList<HspStat> hspStats= new ArrayList<>();
        //replace all hsp in this new one


        for (int i = 0; i<hsps.size(); i++){

            //pour garder que les i index represente la ligne de la database
            if (hsps.get(i).isEmpty()){
                hspStats.add(null);
            }

            for(int j = 0; j< hsps.get(i).size(); j++){



                int startSequenceIndex = sequencePositions.get(i).get(j);
                int startDbIndex = databasePositions.get(i).get(j);

                String currHsp = hsps.get(i).get(j);
                int bruteScore = calculateBruteScore(hsps,database,databasePositions,i,j);
                long bitScore = Math.round((lamda*bruteScore - Math.log(k))/Math.log(2));
                double eValue = databaseTotalSize(database)*sequence.length()/Math.pow(2,bitScore);

                HspStat currHspStat = new HspStat(currHsp,bruteScore,bitScore,eValue,startSequenceIndex,startDbIndex);


                hspStats.add(currHspStat);
            }
        }
        return hspStats;
    }



    public static void bestBitScore (ArrayList<HspStat> hspStats) {

        int maxi = 0;
        for (int i = 0; i < hspStats.size(); i++) {
            for (int j = 0; j < hspStats.size() - 1; j++) {

                if (hspStats.get(i).bitScore > hspStats.get(i + 1).bitScore) {
                    maxi = i;
                }
            }
        }
        System.out.println(hspStats.get(maxi).bitScore);

    }






    //1.7
    public static void main(String[] args) throws FileNotFoundException {
        //System.out.println(readFasta("src/tRNAs.fasta"));
        ArrayList<DBentry> db = new ArrayList<DBentry>();
        ArrayList<DBentry> input = readFasta("src/unknown.fasta");

        //readFasta("src/tRNAs.fasta", database);
        ArrayList<DBentry> unknown = readFasta("src/unknown.fasta");
        ArrayList<DBentry> database = readFasta("src/tRNAs.fasta");
        ArrayList<ArrayList<Integer>> position = new ArrayList<ArrayList<Integer>>();
        ArrayList<ArrayList<Integer>> positionInput = new ArrayList<ArrayList<Integer>>();


        ArrayList<ArrayList<String>> hsps = hsp(unknown.get(0).getSequence(),database , "11111111111", position, positionInput);


        System.out.println("hsp de base");
        System.out.println(hsps);
        //positions
//      System.out.println(position);
//      System.out.println(positionInput);

        System.out.println("extended hsps");
//        ArrayList<ArrayList<String>> extendedHsps = glouton(database,unknown.get(0).getSequence(),positionInput,position,hsps,4,"11111111111");
//        System.out.println(extendedHsps);
        ArrayList<ArrayList<String>> extendedHsps = extendHsps(hsps,database,unknown.get(0).getSequence(),position,positionInput,5);
        System.out.println(extendedHsps);

        //positions
        System.out.println(position);
        System.out.println(positionInput);

        System.out.println("fusion hsp");
        ArrayList<ArrayList<String>> fusion = fusionHsp(extendedHsps,database,position, positionInput);
        System.out.println(fusion);


        //calcul
//        getHspStat(fusion,database,unknown.get(0).getSequence(),position,positionInput);
    }
}
