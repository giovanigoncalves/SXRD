from PIL import Image

images = [
    Image.open("C:/Users/T-Gamer/Downloads/Norma/ASTM_Standard/"
    + f for f in ["1.png", 
                  "2.png",
                  "3.png",
                  "4.png",
                  "5.png",
                  "6.png",
                  "7.png",
                  "8.png",
                  "9.png",
                  "10.png",
                  "11.png",
                  "12.png",
                  "13.png",
                  "14.png"]
            )
    ]

images[0].save("./norma", "PDF", resolution=100.0, save_all=True, append_images=images[1:])