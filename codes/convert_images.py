from sys import argv
from pathlib import Path
import fitz  # PyMuPDF
from PIL import Image

def convert_pdf_to_images(pdf_path: Path, output_folder: Path, output_format='PNG', dpi=300):
    """
    Convert a PDF file to images using PyMuPDF.
    """
    # Create the output folder if it doesn't exist
    output_folder.mkdir(parents=True, exist_ok=True)

    # Open the PDF
    pdf_document = fitz.open(str(pdf_path))

    # Iterate over pages
    for page_num in range(pdf_document.page_count):
        page = pdf_document.load_page(page_num)
        
        # Convert page to a pixmap (image)
        pix = page.get_pixmap(dpi=dpi)

        # Convert pixmap to a PIL Image
        image = Image.frombytes("RGB", [pix.width, pix.height], pix.samples)

        # Save the image
        image_path = output_folder / f"{pdf_path.stem}.{output_format.lower()}"
        image.save(image_path, output_format)
        print(f"Page {page_num + 1} converted and saved as {image_path.name}")

    pdf_document.close()

if __name__ == '__main__':
    # If the given argument is a directory, convert all PDFs in the directory
    # Otherwise, convert the given PDF file
    if Path(argv[1]).is_dir():
        for pdf_path in Path(argv[1]).rglob("*.pdf"):
            convert_pdf_to_images(pdf_path, Path(argv[2]))
    else:
        convert_pdf_to_images(Path(argv[1]), Path(argv[2]))
        
# from sys import argv
# from pathlib import Path
# from pdf2image import convert_from_path
# from PIL import Image

# def convert_pdf_to_images(pdf_path: Path, output_folder: Path, output_format='PNG', dpi=300):
#     """
#     Convert a PDF file to images.
#     """
#     # Create the output folder if it doesn't exist
#     output_folder.mkdir(parents=True, exist_ok=True)

#     # Convert the PDF to images
#     image = convert_from_path(pdf_path, dpi=dpi)


#     # Save the images
#     image_path = output_folder / f"{pdf_path.stem}.{output_format.lower()}"
#     image[0].save(image_path, output_format)
    
#     print(f"PDF converted to image in {image_path.name}")
        
# if __name__ == '__main__':
#     # If the given argument is a directory, convert all PDFs in the directory
#     # Otherwise, convert the given PDF file
#     if Path(argv[1]).is_dir():
#         for pdf_path in Path(argv[1]).rglob("*.pdf"):
#             convert_pdf_to_images(pdf_path, Path(argv[2]))
#     else:
#         convert_pdf_to_images(Path(argv[1]), Path(argv[2]))
        