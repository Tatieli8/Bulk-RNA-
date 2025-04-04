import qupath.lib.gui.scripting.QPEx
import qupath.lib.objects.PathObject
import java.nio.file.Paths
import java.nio.file.Files
import java.nio.charset.StandardCharsets

// Load the current QuPath project
def currentProject = QPEx.getProject()
if (currentProject == null) {
    print("❌ No project loaded!")
    return
}

// Define the output file path
def outputFile = Paths.get("C:/Users/YourUsername/Desktop/Exported_Data.csv")

// Get the list of images in the project
def imageEntries = currentProject.getImageList()
if (imageEntries.isEmpty()) {
    print("❌ No images found in the project!")
    return
}

// Define the CSV headers
def headers = [
    "Image", "Object ID", "Classification", "Centroid X µm", "Centroid Y µm",
    "Nucleus: Area", "Nucleus: Perimeter", "Nucleus: Circularity",
    "Cell: Area", "Cell: Perimeter", "Cell: Circularity"
]

// Initialize list to store data rows
def dataRows = []
dataRows << headers.join(";")

// Loop through all images
for (imgEntry in imageEntries) {
    def imgData = imgEntry.readImageData()
    if (imgData == null) continue // Skip if image cannot be loaded
    def imgHierarchy = imgData.getHierarchy()
    def imageName = imgEntry.getImageName()

    // Get all objects (cells, detections)
    def pathObjects = imgHierarchy.getFlattenedObjectList(null)

    // Collect data for each object
    for (obj in pathObjects) {
        def classification = obj.getPathClass() ? obj.getPathClass().toString() : "Unclassified"
        def centroidX = obj.getROI()?.getCentroidX()  // X coordinate of centroid
        def centroidY = obj.getROI()?.getCentroidY()  // Y coordinate of centroid
        def objId = obj.getID()

        // Extract measurements
        def nucleusArea = obj.getMeasurementList().getMeasurementValue("Nucleus: Area") ?: "N/A"
        def nucleusPerimeter = obj.getMeasurementList().getMeasurementValue("Nucleus: Perimeter") ?: "N/A"
        def nucleusCircularity = obj.getMeasurementList().getMeasurementValue("Nucleus: Circularity") ?: "N/A"
        def cellArea = obj.getMeasurementList().getMeasurementValue("Cell: Area") ?: "N/A"
        def cellPerimeter = obj.getMeasurementList().getMeasurementValue("Cell: Perimeter") ?: "N/A"
        def cellCircularity = obj.getMeasurementList().getMeasurementValue("Cell: Circularity") ?: "N/A"

        // Prepare the data row
        def values = [
            imageName,
            objId,
            classification,
            centroidX != null ? centroidX : "N/A",
            centroidY != null ? centroidY : "N/A",
            nucleusArea,
            nucleusPerimeter,
            nucleusCircularity,
            cellArea,
            cellPerimeter,
            cellCircularity
        ]

        // Add data to the list
        dataRows << values.join(";")
    }
}

// Write data to CSV
Files.write(outputFile, dataRows, StandardCharsets.UTF_8)
print "✅ Data exported successfully to: " + outputFile
