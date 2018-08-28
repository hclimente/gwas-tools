ls in*.gen >gens
ls in*.sample >samples

gtool -M --g gens --s samples --og out.gen --os out.sample --threshold 0.95