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



                int maxScoreLeft = 0;
                int maxScoreRight = 0;
                int currLeftScore = 0;
                int currRightScore = 0;


                boolean leftMatch = false;
                boolean rightMatch = false;



                while ((currLeftScore > maxScoreLeft - seuil || currRightScore > maxScoreRight - seuil) && (!endLeft || !endRight)) {


                    //outOfBoundsVerification
                    if(kmerDBstartPosition-1 == -1 || kmerSequenceStartPosition - 1 == -1){
                        endLeft = true;
                    }
                    if(kmerDBendPosition + 1 == database.get(i).getSequence().length() || kmerSequenceEndPosition + 1  == sequence.length()){
                        endRight = true;
                    }


                    //chek left
                        if (!endLeft && currLeftScore > maxScoreLeft - seuil) {
                            if (database.get(i).getSequence().charAt(kmerDBstartPosition - 1) == sequence.charAt(kmerSequenceStartPosition - 1)) {
                                leftMatch = true;
                                currLeftScore +=5;
                                if (currLeftScore > maxScoreLeft) {
                                    maxScoreLeft = currLeftScore;
                                }
                            } else {
                                currLeftScore -=4;
                            }
                        }


                    //chek right
                    if (!endRight && currRightScore > maxScoreRight - seuil) {
                        if (database.get(i).getSequence().charAt(kmerDBendPosition + 1) == sequence.charAt(kmerSequenceEndPosition + 1)) {
                            rightMatch = true;
                            currRightScore +=5;
                            if(currRightScore>maxScoreRight) {
                                maxScoreRight = currRightScore;
                            }
                        } else {
                            currRightScore -=4;
                        }
                    }


                    //EXTENSION CHOICE -> Based on match -> if both miss -> based on score

                    //extend left (left match && rightMatch or only left match)
                    if (!endLeft && leftMatch) {
                        char leftCar = database.get(i).getSequence().charAt(kmerDBstartPosition - 1);
                        //extend hsp
                        hsps.get(i).set(j, leftCar + hsps.get(i).get(j));
                        //extend start positions
                        dbPositions.get(i).set(j, kmerDBstartPosition - 1);
                        sequencePositions.get(i).set(j, kmerSequenceStartPosition - 1);

                        //next car
                        kmerDBstartPosition--;
                        kmerSequenceStartPosition--;

                    }

                    //extend right (right match only)
                    if (!endRight && rightMatch && !leftMatch){
                        //extend right
                        char rightCar = database.get(i).getSequence().charAt(kmerDBendPosition + 1);
                        //extend hsp
                        hsps.get(i).set(j, hsps.get(i).get(j) + rightCar);


                        //next car
                        kmerDBendPosition++;
                        kmerSequenceEndPosition++;

                    }


                    //if left and right false ->chose better score
                    if (!rightMatch && !leftMatch){


                        if(!endLeft && currLeftScore>= currRightScore){
                            //extend left
                            char leftCar = database.get(i).getSequence().charAt(kmerDBstartPosition - 1);
                            //extend hsp
                            hsps.get(i).set(j, leftCar + hsps.get(i).get(j));
                            //extend start positions
                            dbPositions.get(i).set(j, kmerDBstartPosition - 1);
                            sequencePositions.get(i).set(j, kmerSequenceStartPosition - 1);

                            //next car
                            kmerDBstartPosition--;
                            kmerSequenceStartPosition--;

                        }else if (!endRight){
                            //extend right
                            char rightCar = database.get(i).getSequence().charAt(kmerDBendPosition + 1);
                            //extend hsp
                            hsps.get(i).set(j, hsps.get(i).get(j) + rightCar);


                            //next car
                            kmerDBendPosition++;
                            kmerSequenceEndPosition++;
                        }
                    }

                    leftMatch = false;
                    rightMatch = false;


                }


            }
        }

        return hsps;
    }

    //1.5
    public static ArrayList<ArrayList<String>> fusionHsp (ArrayList<ArrayList<String>> extendedHsps,ArrayList<DBentry> database, ArrayList<ArrayList<Integer>> dbPositions) {

        //condition de fusion : end index dun hsp > start position dun autre && < end
        //jmodifie la hsp sa sera le bordel dedans mais on senfou on ne la plus besoin apres
        //juste besoin de prendre le plusse grosse fusion a la fin

        int startOfMax = 0;
        int endOfMax = 0;
        String fusionMax = "zzz";
        int dbSequence = 0;


            for (int i = 0; i < extendedHsps.size(); i++) {
                for (int j = 0; j < extendedHsps.get(i).size(); j++) {

                    //end index
                    int hspSize = extendedHsps.get(i).get(j).length();
                    int hspEndIndex = dbPositions.get(i).get(j) + hspSize - 1;

                    //chek all other hsps
                    for (int k = 0; k < extendedHsps.size(); k++) {
                        for (int l = 0; k < extendedHsps.get(k).size(); k++) { //tbnk 4 boucle

                            //start and end index of another hsp
                            int otherHspStart = dbPositions.get(k).get(l);
                            int otherHspSize = extendedHsps.get(k).get(l).length();
                            int otherHspEnd = dbPositions.get(k).get(l) + otherHspSize - 1;

                            //condition de fusion
                            if (hspEndIndex > otherHspStart && hspEndIndex < otherHspEnd) {


                                //a partir de end index du current hsp, concat lautre hsp a part sa partie chevauchee
                                int longueurChevauchement = hspEndIndex - otherHspStart;

                                String concatPartOfOtherHsp = extendedHsps.get(k).get(l).substring(longueurChevauchement);
                                String fusionFinale = extendedHsps.get(i).get(j) + concatPartOfOtherHsp;

                                //get les index de debut et de fin de la sequence fusionee
                                int start = dbPositions.get(i).get(j);
                                int end = start + fusionFinale.length()-1;
                                //add fusion
                                extendedHsps.get(i).set(j,start + " " +fusionFinale + " " + end);


                                //Fusion Max ?
                                if (fusionFinale.length()> fusionMax.length()){
                                    fusionMax = fusionFinale;
                                    startOfMax = dbPositions.get(i).get(j);
                                    endOfMax = startOfMax + fusionFinale.length()-1;

                                    //database line
                                    dbSequence = i;
//                                    System.out.println(i);
                                }

                            }

                        }

                    }

                }
            }

            System.out.println(startOfMax + " " + fusionMax + " " + endOfMax);

        //}

        return extendedHsps;
    }








        //1.6
    //1.7
    public static void main(String[] args) throws FileNotFoundException {
        //System.out.println(readFasta("src/tRNAs.fasta"));
        ArrayList<DBentry> db = new ArrayList<DBentry>();
        ArrayList<DBentry> input = readFasta("src/unknown.fasta");
        DBentry x = new DBentry("BRUH","AAAB");
        db.add(x);
        DBentry y = new DBentry("BRUH1","BAAA");
        db.add(y);
        DBentry z = new DBentry("BRUH2","ACAA");
        db.add(z);
        //readFasta("src/tRNAs.fasta", database);
        ArrayList<DBentry> unknown = readFasta("src/unknown.fasta");
        ArrayList<DBentry> database = readFasta("src/tRNAs.fasta");
        ArrayList<ArrayList<Integer>> position = new ArrayList<ArrayList<Integer>>();
        ArrayList<ArrayList<Integer>> positionInput = new ArrayList<ArrayList<Integer>>();
        //System.out.println(kmer(unknown.get(1), 11));

        //Creer ici la database pour avoir acces partout
        //System.out.println("GGGAGAATGACTGAGTGGTTAAAAGTGACAGACTGTAAATCTGTTGAAATTATTTCTACGTAGGTTCGAATCCTGCTTCTCCCA".length());

        //HSP
        //System.out.println(hsp(unknown.get(0).getSequence(),database , "11111111111", position, positionInput));

        ArrayList<ArrayList<String>> hsps = hsp(unknown.get(0).getSequence(),database , "11111111111", position, positionInput);

        System.out.println("hsp de base");
        System.out.println(hsps);
        //positions
//      System.out.println(position);
//      System.out.println(positionInput);

        System.out.println("extended hsps");
        ArrayList<ArrayList<String>> extendedHsps = extendHsps(hsps,database,unknown.get(0).getSequence(),position,positionInput,5);
        System.out.println(extendedHsps);

        //positions
        System.out.println(position);
        System.out.println(positionInput);

        System.out.println("fusion hsp");
        System.out.println(fusionHsp(extendedHsps,database,position));


    }
}
