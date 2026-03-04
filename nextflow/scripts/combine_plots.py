#!/usr/bin/env python3
"""Combine two Manhattan plots vertically"""

import sys
from PIL import Image, ImageDraw, ImageFont

def combine_images_vertical(img1_path, img2_path, output_path, label1="A", label2="B"):
    """Combine two images vertically with labels"""
    
    print("="*60)
    print("Combining Manhattan Plots Vertically")
    print("="*60)
    print(f"Top panel:    {img1_path}")
    print(f"Bottom panel: {img2_path}")
    print(f"Output:       {output_path}")
    print(f"Labels:       {label1} | {label2}\n")
    
    # Load images
    print("Reading images...")
    img1 = Image.open(img1_path)
    img2 = Image.open(img2_path)
    
    print(f"  Plot 1: {img1.width} x {img1.height} pixels")
    print(f"  Plot 2: {img2.width} x {img2.height} pixels\n")
    
    # Make sure both have the same width
    max_width = max(img1.width, img2.width)
    if img1.width != max_width:
        img1 = img1.resize((max_width, int(img1.height * max_width / img1.width)), Image.LANCZOS)
    if img2.width != max_width:
        img2 = img2.resize((max_width, int(img2.height * max_width / img2.width)), Image.LANCZOS)
    
    # Create combined image
    print("Combining images...")
    total_height = img1.height + img2.height
    combined = Image.new('RGB', (max_width, total_height))
    
    # Paste images
    combined.paste(img1, (0, 0))
    combined.paste(img2, (0, img1.height))
    
    # Add labels
    print("Adding labels...")
    draw = ImageDraw.Draw(combined)
    try:
        font = ImageFont.truetype("/usr/share/fonts/dejavu/DejaVuSans-Bold.ttf", 60)
    except:
        font = ImageFont.load_default()
    
    # Add label to top panel
    draw.text((40, 40), label1, fill='black', font=font)
    # Add label to bottom panel  
    draw.text((40, img1.height + 40), label2, fill='black', font=font)
    
    # Save
    print("Saving combined image...")
    combined.save(output_path, 'PNG', quality=95)
    
    file_size = os.path.getsize(output_path) / (1024**2)
    print(f"\n{'='*60}")
    print(f"Combined plot saved: {output_path}")
    print(f"Final dimensions: {combined.width} x {combined.height} pixels")
    print(f"File size: {file_size:.2f} MB")
    print("="*60)

if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("Usage: python combine_plots.py <plot1.png> <plot2.png> <output.png> [label1] [label2]")
        sys.exit(1)
    
    import os
    
    plot1 = sys.argv[1]
    plot2 = sys.argv[2]
    output = sys.argv[3]
    label1 = sys.argv[4] if len(sys.argv) > 4 else "A"
    label2 = sys.argv[5] if len(sys.argv) > 5 else "B"
    
    combine_images_vertical(plot1, plot2, output, label1, label2)
