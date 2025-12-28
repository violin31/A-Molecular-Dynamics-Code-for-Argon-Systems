if [ "$#" -ne 3 ]
then
    # Not enough input arguments
    echo "Usage compile_and_run.sh <file_path> <delta_time> <max_steps>"
    echo "<file_path> The input file containing the data for the molecule"
    echo "<delta_time> Parameter of the Verlet Algorthm"
    echo "<max_steps> Number of iterations of the Verlet Algorithm"
    echo "Exit"
    exit
fi

file_path=$1
delta_time=$2
max_steps=$3

mkdir -p lib
gcc -I/usr/local/include src/dynamics.c -lm -o ./lib/main 
export LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH
./lib/main ${file_path} ${delta_time} ${max_steps}