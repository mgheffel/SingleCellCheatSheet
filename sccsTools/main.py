import argparse
import scanpy as sc

def main(adata_path, umap_color):
    # Read the AnnData object
    adata = sc.read_h5ad(adata_path)

    # Create a UMAPPlot object
    umap_plot = sc.pl.umap(adata,color=umap_color,save='out.png')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate UMAP plot from AnnData object and save as PNG.')
    parser.add_argument('adata_path', type=str, help='Path to the AnnData file (.h5ad).')
    parser.add_argument('umap_color', type=str, help='Color to use for the UMAP plot.')

    args = parser.parse_args()
    main(args.adata_path, args.umap_color)
