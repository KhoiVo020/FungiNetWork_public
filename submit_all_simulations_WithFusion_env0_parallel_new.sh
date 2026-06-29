#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --cpus-per-task=1
#SBATCH --mem=200G
#SBATCH --time=480:00:00
#SBATCH --mail-user=kvo020@ucr.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name=WithFusion_env3_parallel
#SBATCH --output=/bigdata/alberlab/kvo020/batch_WithFusion_env3_%j.out
#SBATCH --error=/bigdata/alberlab/kvo020/batch_WithFusion_env3_%j.err
#SBATCH -p intel

module load anaconda
conda activate Khoi_env

DIRECTORIES=(
    "/bigdata/alberlab/kvo020/Environment0/Hyphae_WithFusion_01"
    "/bigdata/alberlab/kvo020/Environment0/Hyphae_WithFusion_02"
    "/bigdata/alberlab/kvo020/Environment0/Hyphae_WithFusion_03"
    "/bigdata/alberlab/kvo020/Environment0/Hyphae_WithFusion_04"
    "/bigdata/alberlab/kvo020/Environment0/Hyphae_WithFusion_05"
    "/bigdata/alberlab/kvo020/Environment0/Hyphae_WithFusion_06"
    "/bigdata/alberlab/kvo020/Environment0/Hyphae_WithFusion_07"
    "/bigdata/alberlab/kvo020/Environment0/Hyphae_WithFusion_08"
    "/bigdata/alberlab/kvo020/Environment0/Hyphae_WithFusion_09"
    "/bigdata/alberlab/kvo020/Environment0/Hyphae_WithFusion_10"
)
s
# Update cwd_path in each driver file before running
for DIR in "${DIRECTORIES[@]}"; do
    DRIVER_FILE="${DIR}/driver_fungalGrowth_singleNutrient.py"
    if [ -f "$DRIVER_FILE" ]; then
        sed -i "s|^cwd_path.*=.*|cwd_path = '${DIR}'|" "$DRIVER_FILE"
        echo "Updated cwd_path in $DRIVER_FILE"
    else
        echo "Warning: $DRIVER_FILE not found"
    fi
done

# Run all simulations in parallel
for DIR in "${DIRECTORIES[@]}"; do
    echo "Starting simulation in: $DIR"
    cd "$DIR"
    python3 "${DIR}/driver_fungalGrowth_singleNutrient.py" -i "${DIR}/parameters.ini" &
done

# Wait for all background jobs to complete
wait

echo "All simulations completed!"
conda deactivate
