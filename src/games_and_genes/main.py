from games_and_genes.configs import cfg

def run():
    print(f"Running experiment: {cfg['name']}")
    print(f"Data directory: {cfg['data_dir']}")
    print(f"Output directory: {cfg['output_dir']}")

if __name__ == "__main__":
    run()